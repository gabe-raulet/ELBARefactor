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

bool KmerParsePass
(
    Set<TKmer>&            localkmers, /* received k-mers will go here */
    const Vector<String>&  myreads,    /* local read sequence set */
    size_t&                myoffset,   /* current read sequence offset */
    Bloom&                 bm,         /* Bloom filter */
    SharedPtr<CommGrid>    commgrid    /* communicator info */
)
{
    /*
     * Each process participates each round and needs to know itself.
     */
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    int len;
    size_t nreads = myreads.size(); /* number of locally stored reads */
    size_t maxsending = 0;
    size_t kmersthisbatch = 0;

    size_t bytesperkmer = TKmer::N_BYTES; /* generally 16 bytes needed per k-mer because 32<k<64 means 64-128 bits, so 128 bits necessary = 2^7/2^3 bytes = 16 bytes */
    size_t bytesperentry = bytesperkmer + 4; /* each entry needs 4 more bytes (apparently) */
    size_t memthreshold = (MAX_ALLTOALL_MEM / nprocs) * 2; /* set memory threshold to 2*(128 MB / nprocs) (?) */

    Vector<TKmer> mykmers; /* vector of k-mers I will receive (should probably be declared much later in this method) */
    Vector<Vector<TKmer>> outgoing(nprocs); /* a bucket for each processor destination to put k-mers in */

    /*
     * Remember this function is called repeatedly
     * in a loop at the higher stack frame. The myoffset variable
     * is basically a local variable in that higher stack
     * frame and thus a 'global variable' within this function (KmerParsePass),
     * so it is constantly moving towards nreads. The reason for
     * the strange semantics here is because we have to repeatedly
     * stage communication phases bounded by memory, which requires
     * repeatedly recomputing new parameters for the communication
     * functions each phase.
     */
    while (myoffset < nreads)
    {
        /*
         * We get the next read stored locally;
         */
        String myseq = myreads[myoffset++];
        len = myseq.size();

        if (len <= TKmer::k) /* expect read to be longer than k */
            continue;

        Vector<TKmer> kmers = TKmer::GetRepKmers(myseq); /* get every representative k-mer in order from read */
        kmersthisbatch += kmers.size(); /* we're sending every k-mer so we increment the kmersthisbatch parameter */

        /* go through each k-mer */
        for (size_t j = 0; j < kmers.size(); ++j)
        {
            int owner = kmers[j].GetInferredOwner(nprocs); /* compute processor destination id for kmers[j] */
            outgoing[owner].push_back(kmers[j]); /* add the k-mer to the outgoing bucket corresponding to proccessor destination id */
            size_t sending = outgoing[owner].size(); /* get size of most recently augmented bucket */
            maxsending = std::max(sending, maxsending); /* update maxsending to be the largest bucket seen */
        }

        /*
         * if we don't have enough memory left ?? (I'm starting to think the memory limitation here
         * is nothing more than the MPI-3 2^32-1 counts limitation, so perhaps we can get rid of
         * these repeated calls using MPI-4 2^64-1 counts) then break out of loop to go to
         * communicaton parameter computation (and then communictaion)
         */
        if (maxsending * bytesperentry >= memthreshold || (kmersthisbatch + len) * bytesperentry >= MAX_ALLTOALL_MEM)
            break;
    }

    bool finished = (myoffset >= nreads); /* this processor will be finished reading its local reads if this condition is met */

    /*
     * Alltoallv communication parameters.
     */
    std::vector<int> sendcnt(nprocs);
    std::vector<int> recvcnt(nprocs);
    std::vector<int> sdispls(nprocs);
    std::vector<int> rdispls(nprocs);

    /*
     * Each processor is going to be sent a specific number of k-mers.
     */
    for (int i = 0; i < nprocs; ++i)
    {
        sendcnt[i] = outgoing[i].size() * bytesperkmer;
    }

    /*
     * Share the communication paramters.
     */
    MPI_Alltoall(sendcnt.data(), 1, MPI_INT, recvcnt.data(), 1, MPI_INT, commgrid->GetWorld());

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
     * setup send buffer for filling.
     */
    grow_buffer(scratch1, totsend);

    uint8_t *sendbuf = (uint8_t*) get_start_buffer(scratch1);

    /*
     * for each processor i
     */
    for (int i = 0; i < nprocs; ++i)
    {
        size_t nkmers2send = sendcnt[i]; /* number of k-mers going to proc i */
        uint8_t *addrs2fill = sendbuf + sdispls[i]; /* get location within buffer to fill for proc i */

        /* for each k-mer */
        for (size_t j = 0; j < nkmers2send; ++j)
        {
            (outgoing[i][j]).CopyDataInto(addrs2fill); /* */
            addrs2fill += bytesperkmer;
        }

        outgoing[i].clear();
    }

    grow_buffer(scratch2, totsend);

    uint8_t *recvbuf = (uint8_t*)get_start_buffer(scratch2);

    MPI_Alltoallv(sendbuf, sendcnt.data(), sdispls.data(), MPI_BYTE, recvbuf, recvcnt.data(), rdispls.data(), MPI_BYTE, commgrid->GetWorld());

    uint64_t nkmersrecvd = totrecv / bytesperkmer;

    for (uint64_t i = 0; i < nkmersrecvd; ++i)
    {
        TKmer kk;
        kk.CopyDataFrom(recvbuf + (i * bytesperkmer));
        mykmers.push_back(kk);
    }

    for (auto itr = mykmers.begin(); itr != mykmers.end(); ++itr)
    {
        if (bm.Check(itr->GetBytes(), TKmer::N_BYTES))
        {
            if (localkmers.find(*itr) == localkmers.end())
                localkmers.insert(*itr);
        }
        else
        {
            bm.Add(itr->GetBytes(), TKmer::N_BYTES);
        }
    }

    mykmers.clear();

    int isfinished = static_cast<int>(finished);

    MPI_Allreduce(MPI_IN_PLACE, &isfinished, 1, MPI_INT, MPI_LAND, commgrid->GetWorld());

    return static_cast<bool>(isfinished);
}

Set<TKmer> GetLocalKmers(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid)
{
    Set<TKmer> localkmers;

    /*
     * Each process participates each round and needs to know itself.
     */
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    double cardinality = EstimateKmerCardinality(myreads, commgrid);

    if (!myrank) std::cerr << "Estimated k-mer  cardinality is " << cardinality << std::endl;

    size_t amt = static_cast<size_t>(cardinality * 0.1 + 64000); /* ?? */
    int64_t intcard = static_cast<int64_t>(std::ceil(cardinality));

    Bloom bm(intcard, 0.05);

    scratch1 = init_buffer(MAX_ALLTOALL_MEM);
    scratch2 = init_buffer(MAX_ALLTOALL_MEM);

    size_t myoffset = 0;
    size_t numreads = myreads.size();

    bool finished;

    int len;
    size_t nreads = myreads.size(); /* number of locally stored reads */
    size_t maxsending = 0;
    size_t kmersthisbatch = 0;

    size_t bytesperkmer = TKmer::N_BYTES; /* generally 16 bytes needed per k-mer because 32<k<64 means 64-128 bits, so 128 bits necessary = 2^7/2^3 bytes = 16 bytes */
    size_t bytesperentry = bytesperkmer + 4; /* each entry needs 4 more bytes (apparently) */

    Vector<TKmer> mykmers; /* vector of k-mers I will receive (should probably be declared much later in this method) */
    Vector<Vector<TKmer>> outgoing(nprocs); /* a bucket for each processor destination to put k-mers in */

    while (myoffset < nreads)
    {
        /*
         * We get the next read stored locally;
         */
        String myseq = myreads[myoffset++];
        len = myseq.size();

        if (len <= TKmer::k) /* expect read to be longer than k */
            continue;

        Vector<TKmer> kmers = TKmer::GetRepKmers(myseq); /* get every representative k-mer in order from read */
        kmersthisbatch += kmers.size(); /* we're sending every k-mer so we increment the kmersthisbatch parameter */

        /* go through each k-mer */
        for (size_t j = 0; j < kmers.size(); ++j)
        {
            int owner = kmers[j].GetInferredOwner(nprocs); /* compute processor destination id for kmers[j] */
            outgoing[owner].push_back(kmers[j]); /* add the k-mer to the outgoing bucket corresponding to proccessor destination id */
            size_t sending = outgoing[owner].size(); /* get size of most recently augmented bucket */
            maxsending = std::max(sending, maxsending); /* update maxsending to be the largest bucket seen */
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
       sendcnt[i] = outgoing[i].size() * bytesperkmer;
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
    * setup send buffer for filling.
    */
   grow_buffer(scratch1, totsend);

   uint8_t *sendbuf = (uint8_t*) get_start_buffer(scratch1);

   /*
    * for each processor i
    */
   for (int i = 0; i < nprocs; ++i)
   {
       size_t nkmers2send = sendcnt[i]; /* number of k-mers going to proc i */
       uint8_t *addrs2fill = sendbuf + sdispls[i]; /* get location within buffer to fill for proc i */

       /* for each k-mer */
       for (size_t j = 0; j < nkmers2send; ++j)
       {
           (outgoing[i][j]).CopyDataInto(addrs2fill); /* copy next k-mer from bucket into buffer (addres2fill) */
           addrs2fill += bytesperkmer; /* increment buffer pointer */
       }

       outgoing[i].clear(); /* don't need it anymore */
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
