#include "KmerOps.h"
#include "Bloom.h"
#include "Buffer.h"
#include <numeric>
#include <algorithm>
#include <cmath>

buffer_t *scratch1 = NULL;
buffer_t *scratch2 = NULL;

double EstimateKmerCardinality(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid)
{
    HyperLogLog hll(12);

    for (auto itr = myreads.begin(); itr != myreads.end(); ++itr)
    {
        TKmer::InsertIntoHLL(*itr, hll);
    }

    hll.ParallelMerge(commgrid->GetWorld());

    return hll.Estimate();
}

Set<TKmer> GetLocalKmers(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    Set<TKmer> localkmers;

    double cardinality = EstimateKmerCardinality(myreads, commgrid);

    size_t amt = static_cast<size_t>(cardinality * 0.1 + 64000); /* ?? */
    int64_t intcard = static_cast<int64_t>(std::ceil(cardinality));

    Bloom bm(intcard, 0.05);

    scratch1 = init_buffer(MAX_ALLTOALL_MEM);
    scratch2 = init_buffer(MAX_ALLTOALL_MEM);

    size_t numreads = myreads.size();

    size_t bytesperkmer = TKmer::N_BYTES; /* generally 16 bytes needed per k-mer because 32<k<64 means 64-128 bits, so 128 bits necessary = 2^7/2^3 bytes = 16 bytes */
    size_t bytesperentry = bytesperkmer + 4; /* each entry needs 4 more bytes (apparently) */

    Vector<TKmer> mykmers; /* vector of k-mers I will receive (should probably be declared much later in this method) */
    Vector<Vector<TKmer>> outgoing_kmers(nprocs); /* a bucket for each processor destination to put k-mers in */

    for (auto readitr = myreads.begin(); readitr != myreads.end(); ++readitr)
    {
        if (readitr->size() <= TKmer::k) /* expect read to be longer than k */
            continue;

        Vector<TKmer> kmers = TKmer::GetRepKmers(*readitr); /* get every representative k-mer in order from read */

        for (auto meritr = kmers.begin(); meritr != kmers.end(); ++meritr)
        {
            int owner = meritr->GetInferredOwner(nprocs); /* compute processor destination id for kmers[j] */
            outgoing_kmers[owner].push_back(*meritr); /* add the k-mer to the outgoing bucket corresponding to proccessor destination id */
        }
    }

   /*
    * Alltoallv communication parameters.
    */
   std::vector<MPI_Count> sendcnt(nprocs);
   std::vector<MPI_Count> recvcnt(nprocs);
   std::vector<MPI_Aint> sdispls(nprocs);
   std::vector<MPI_Aint> rdispls(nprocs);

   /*
    * Each processor is going to be sent a specific number of k-mers.
    */
   for (int i = 0; i < nprocs; ++i)
   {
       sendcnt[i] = outgoing_kmers[i].size() * bytesperkmer;
   }

   /*
    * Share the communication paramters.
    */
   MPI_Alltoall(sendcnt.data(), 1, MPI_COUNT, recvcnt.data(), 1, MPI_COUNT, commgrid->GetWorld());

   /*
    * Comupte displacement vectors.
    */
   sdispls.front() = 0;
   rdispls.front() = 0;

   std::partial_sum(sendcnt.begin(), sendcnt.end()-1, sdispls.begin()+1);
   std::partial_sum(recvcnt.begin(), recvcnt.end()-1, rdispls.begin()+1);

   int64_t totsend = std::accumulate(sendcnt.begin(), sendcnt.end(), 0);
   int64_t totrecv = std::accumulate(recvcnt.begin(), recvcnt.end(), 0);

   /*
    * Now that we know exactly how many k-mers we are sending,
    * we can reserve the amount of memory required for a our
    * send buffer and then fill it up.
    */

   grow_buffer(scratch1, totsend);

   uint8_t *sendbuf = (uint8_t*) get_start_buffer(scratch1);

   for (int i = 0; i < nprocs; ++i)
   {
       size_t nkmers2send = sendcnt[i]; /* number of k-mers going to proc i */
       uint8_t *addrs2fill = sendbuf + sdispls[i]; /* get location within buffer to fill for proc i */

       for (size_t j = 0; j < nkmers2send; ++j)
       {
           (outgoing_kmers[i][j]).CopyDataInto(addrs2fill); /* copy next k-mer from bucket into buffer (addres2fill) */
           addrs2fill += bytesperkmer; /* increment buffer pointer */
       }

       outgoing_kmers[i].clear(); /* don't need it anymore */
   }

   grow_buffer(scratch2, totsend);

   uint8_t *recvbuf = (uint8_t*)get_start_buffer(scratch2);

   /*
    * Perform alltoallv communication
    */
   MPI_Alltoallv_c(sendbuf, sendcnt.data(), sdispls.data(), MPI_BYTE, recvbuf, recvcnt.data(), rdispls.data(), MPI_BYTE, commgrid->GetWorld());

   /* number of k-mers received is tot number of bytes received / bytes per k-mer */
   uint64_t nkmersrecvd = totrecv / bytesperkmer;
   uint64_t *dst = (uint64_t *)recvbuf;

   /* reconstruct all the received k-mers and put them on mykmers */
   for (uint64_t i = 0; i < nkmersrecvd; ++i)
   {
       TKmer kk;
       kk.CopyDataFrom(recvbuf + (i * bytesperkmer));
       mykmers.push_back(kk);
       dst += bytesperkmer;
   }

   /* go through each received k-mer */
   for (auto itr = mykmers.begin(); itr != mykmers.end(); ++itr)
   {
       if (bm.Check(itr->GetBytes(), TKmer::N_BYTES))
       {
           /*
            * k-mer was in the bloom filter, which
            * means it has probably been seen before and
            * is therefore probably not unique, so we
            * add it to the hash table if it isn't there
            * (checking the hash table is much more expensive
            * then the Bloom filter)
            */
           if (localkmers.find(*itr) == localkmers.end())
               localkmers.insert(*itr);
       }
       else
       {
           /*
            * k-mer wasn't in the Bloom filter, which
            * means it definitley hasn't been seen yet. Add
            * it to the Bloom filter, so that if it never
            * shows up again we don't have to check the
            * more expensive hash table for a unique k-mer.
            * Remember that a signifncat fraction of k-mers
            * are unique (due to read errors) so this will
            * save significant time during the k-mer discovery
            * phase.
            */
           bm.Add(itr->GetBytes(), TKmer::N_BYTES);
       }
   }

    free_buffer(scratch1);
    free_buffer(scratch2);

    return localkmers;
}

template <bool COUNT_PHASE>
void KmerPass(const Vector<String>& myreads, Map<TKmer, KmerCountType>& kmercounts, SharedPtr<CommGrid> commgrid)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    double cardinality = 0;
    size_t amt = 0;
    int64_t intcard = 0;

    Bloom *bm;

    if (!COUNT_PHASE)
    {
        cardinality = EstimateKmerCardinality(myreads, commgrid);
        amt = static_cast<size_t>(cardinality * 0.1 + 64000); /* ?? */
        intcard = static_cast<int64_t>(std::ceil(cardinality));
        bm = new Bloom(intcard, 0.05);
    }

    scratch1 = init_buffer(MAX_ALLTOALL_MEM);
    scratch2 = init_buffer(MAX_ALLTOALL_MEM);

    size_t numreads = myreads.size();

    size_t bytesperkmer = TKmer::N_BYTES; /* generally 16 bytes needed per k-mer because 32<k<64 means 64-128 bits, so 128 bits necessary = 2^7/2^3 bytes = 16 bytes */
    size_t bytesperentry = bytesperkmer + 4; /* each entry needs 4 more bytes (apparently) */

    Vector<TKmer> mykmers; /* vector of k-mers I will receive (should probably be declared much later in this method) */
    Vector<Vector<TKmer>> outgoing_kmers(nprocs); /* a bucket for each processor destination to put k-mers in */

    Vector<PosInRead> mypositions;
    Vector<ReadId> myreadids;
    Vector<Vector<PosInRead>> outgoing_positions;
    Vector<Vector<ReadId>> outgoing_readids;

    if (COUNT_PHASE)
    {
        outgoing_positions.resize(nprocs);
        outgoing_readids.resize(nprocs);
    }

    ReadId readid = 1;

    for (auto readitr = myreads.begin(); readitr != myreads.end(); ++readitr)
    {
        if (readitr->size() <= TKmer::k) /* expect read to be longer than k */
            continue;

        Vector<TKmer> kmers = TKmer::GetRepKmers(*readitr); /* get every representative k-mer in order from read */

        PosInRead pos = 0;

        for (auto meritr = kmers.begin(); meritr != kmers.end(); ++meritr)
        {
            int owner = meritr->GetInferredOwner(nprocs); /* compute processor destination id for kmers[j] */
            outgoing_kmers[owner].push_back(*meritr); /* add the k-mer to the outgoing bucket corresponding to proccessor destination id */

            if (COUNT_PHASE)
            {
                outgoing_positions[owner].push_back(pos);
                outgoing_readids[owner].push_back(readid);
            }

            pos++;
        }

        readid++;
    }

    /*
     * Alltoallv communication parameters.
     */
    std::vector<MPI_Count> sendcnt(nprocs);
    std::vector<MPI_Count> recvcnt(nprocs);
    std::vector<MPI_Aint> sdispls(nprocs);
    std::vector<MPI_Aint> rdispls(nprocs);

    bytesperentry = bytesperkmer + (COUNT_PHASE? sizeof(ReadId) + sizeof(PosInRead) : 0);

    /*
     * Each processor is going to be sent a specific number of k-mers.
     */
    for (int i = 0; i < nprocs; ++i)
    {
        sendcnt[i] = outgoing_kmers[i].size() * bytesperentry;
    }

    /*
     * Share the communication paramters.
     */
    MPI_Alltoall(sendcnt.data(), 1, MPI_COUNT, recvcnt.data(), 1, MPI_COUNT, commgrid->GetWorld());

    /*
     * Comupte displacement vectors.
     */
    sdispls.front() = 0;
    rdispls.front() = 0;

    std::partial_sum(sendcnt.begin(), sendcnt.end()-1, sdispls.begin()+1);
    std::partial_sum(recvcnt.begin(), recvcnt.end()-1, rdispls.begin()+1);

    int64_t totsend = std::accumulate(sendcnt.begin(), sendcnt.end(), 0);
    int64_t totrecv = std::accumulate(recvcnt.begin(), recvcnt.end(), 0);

    /*
     * Now that we know exactly how many k-mers we are sending,
     * we can reserve the amount of memory required for a our
     * send buffer and then fill it up.
     */

    grow_buffer(scratch1, totsend);

    uint8_t *sendbuf = (uint8_t*) get_start_buffer(scratch1);

    for (int i = 0; i < nprocs; ++i)
    {
        size_t nkmers2send = sendcnt[i]; /* number of k-mers going to proc i */
        uint8_t *addrs2fill = sendbuf + sdispls[i]; /* get location within buffer to fill for proc i */

        for (size_t j = 0; j < nkmers2send; ++j)
        {
            (outgoing_kmers[i][j]).CopyDataInto(addrs2fill); /* copy next k-mer from bucket into buffer (addres2fill) */

            if (COUNT_PHASE)
            {
                *((ReadId *)(addrs2fill + bytesperkmer)) = outgoing_readids[i][j];
                *((PosInRead *)(addrs2fill + bytesperkmer + sizeof(ReadId))) = outgoing_positions[i][j];
            }

            addrs2fill += bytesperentry; /* increment buffer pointer */
        }

        outgoing_kmers[i].clear(); /* don't need it anymore */

        if (COUNT_PHASE)
        {
            outgoing_readids[i].clear();
            outgoing_positions[i].clear();
        }
    }


    grow_buffer(scratch2, totsend);

    uint8_t *recvbuf = (uint8_t*)get_start_buffer(scratch2);

    /*
     * Perform alltoallv communication
     */
    MPI_Alltoallv_c(sendbuf, sendcnt.data(), sdispls.data(), MPI_BYTE, recvbuf, recvcnt.data(), rdispls.data(), MPI_BYTE, commgrid->GetWorld());

    /* number of k-mers received is tot number of bytes received / bytes per k-mer */
    uint64_t nkmersrecvd = totrecv / bytesperentry;
    uint8_t *addrs2read = recvbuf;

    /* reconstruct all the received k-mers and put them on mykmers */
    for (uint64_t i = 0; i < nkmersrecvd; ++i)
    {
        TKmer kk;
        kk.CopyDataFrom(addrs2read);
        mykmers.push_back(kk);

        if (COUNT_PHASE)
        {
            ReadId readid = *((ReadId *)(addrs2read + bytesperkmer));
            PosInRead pos = *((PosInRead *)(addrs2read + bytesperkmer + sizeof(ReadId)));
            myreadids.push_back(readid);
            mypositions.push_back(pos);
        }

        addrs2read += bytesperentry;
    }

    size_t id = 0;

    /* go through each received k-mer */
    for (auto itr = mykmers.begin(); itr != mykmers.end(); ++itr, ++id)
    {
        if (!COUNT_PHASE)
        {
            if (bm->Check(itr->GetBytes(), TKmer::N_BYTES))
            {
                /*
                 * k-mer was in the bloom filter, which
                 * means it has probably been seen before and
                 * is therefore probably not unique, so we
                 * add it to the hash table if it isn't there
                 * (checking the hash table is much more expensive
                 * then the Bloom filter)
                 */
                if (kmercounts.find(*itr) == kmercounts.end())
                {
                    kmercounts.insert({*itr, std::make_tuple(READIDS{}, POSITIONS{}, 0)});
                }
            }
            else
            {
                /*
                 * k-mer wasn't in the Bloom filter, which
                 * means it definitley hasn't been seen yet. Add
                 * it to the Bloom filter, so that if it never
                 * shows up again we don't have to check the
                 * more expensive hash table for a unique k-mer.
                 * Remember that a signifncat fraction of k-mers
                 * are unique (due to read errors) so this will
                 * save significant time during the k-mer discovery
                 * phase.
                 */
                bm->Add(itr->GetBytes(), TKmer::N_BYTES);
            }
        }
        else
        {
            auto got = kmercounts.find(*itr);

            /* if k-mer was parsed in first pass */
            if (got != kmercounts.end())
            {
                READIDS& readids = std::get<0>(got->second);
                POSITIONS& positions = std::get<1>(got->second);
                int& cnt = std::get<2>(got->second);

                if (cnt == UPPER_KMER_FREQ)
                {
                    kmercounts.erase(got);
                }
                else
                {
                    for (int j = 0; j < UPPER_KMER_FREQ; ++j)
                    {
                        if (readids[j] == 0)
                        {
                            readids[j] = myreadids[id];
                            positions[j] = mypositions[id];
                            ++cnt;
                            ++id;
                            break;
                        }
                    }
                }
            }
        }
    }

    free_buffer(scratch1);
    free_buffer(scratch2);

    if (!COUNT_PHASE) delete bm;
}

template void KmerPass<false>(const Vector<String>& myreads, Map<TKmer, KmerCountType>& kmercounts, SharedPtr<CommGrid> commgrid);
template void KmerPass<true>(const Vector<String>& myreads, Map<TKmer, KmerCountType>& kmercounts, SharedPtr<CommGrid> commgrid);

