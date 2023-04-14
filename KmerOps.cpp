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

Set<TKmer> GetKmersSmart(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid)
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
