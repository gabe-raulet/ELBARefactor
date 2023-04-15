#include "KmerComm.h"
#include <numeric>
#include <algorithm>

KmerCountMap GetKmerCountMapKeys(const Vector <String>& myreads, SharedPtr<CommGrid> commgrid)
{
    KmerCountMap kmermap;

    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    size_t numreads = myreads.size();

    Vector<Vector<TKmer>> kmerbuckets(nprocs); /* outgoing k-mer buckets */

    /*
     * When we use the word "k-mer" there are often actually two different things
     * we are referring to, which can lead to confusion. One thing we're interested
     * in is simply that a k-mer exists in a certain read at a certain position.
     * Let's call them "seed k-mers". The other kind of k-mer is simply a specific
     * sequence of k nucleotides. This latter kind of k-mer can appear multiple times
     * in different reads and even in the same read. We still refer to it as one k-mer,
     * because in this case we only care about its particular sequence. We will
     * call these just "k-mers". Hence, multiple distinct "seed k-mers" can in actual
     * fact be the same "k-mer".
     */

    /*
     * Go through each local read.
     */
    for (auto readitr = myreads.begin(); readitr != myreads.end(); ++readitr)
    {
        /*
         * If it is too small then continue to the next one.
         */
        if (readitr->size() < TKmer::k)
            continue;

        /*
         * Get all the representative k-mer seeds.
         */
        Vector<TKmer> repmers = TKmer::GetRepKmers(*readitr);

        /*
         * Go through each k-mer seed.
         */
        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr)
        {
            /*
             * Get destination processor rank of this k-mer seed
             * and put it in the corresponding bucket.
             */
            kmerbuckets[meritr->GetInferredOwner(nprocs)].push_back(*meritr);
        }
    }

    /*
     * Alltoallv communication of k-mers requires communication parameters
     * sendcnt, recvcnt, sdispls, and rdispls.
     */
    std::vector<MPI_Count_type> sendcnt(nprocs); /* sendcnt[i] is number of bytes of k-mers this process sends to process i */
    std::vector<MPI_Count_type> recvcnt(nprocs); /* recvnct[i] is number of bytes of k-mers this process receives from process i */
    std::vector<MPI_Displ_type> sdispls(nprocs); /* sdispls[i] = sdispls[i-1] + sendcnt[i] */
    std::vector<MPI_Displ_type> rdispls(nprocs); /* rdispls[i] = rdispls[i-1] + recvcnt[i] */

    /*
     * Initialize sendcnt parameter for local process.
     */
    for (int i = 0; i < nprocs; ++i)
    {
        /*
         * Each k-mer is a fixed number of bytes (TKmer::N_BYTES),
         * usually 16 for 32 < k < 64.
         */
        sendcnt[i] = kmerbuckets[i].size() * TKmer::N_BYTES;
    }

    /*
     * Let every process know how many k-mers it will be receiving
     * from each other process.
     */
    MPI_ALLTOALL(sendcnt.data(), 1, MPI_COUNT_TYPE, recvcnt.data(), 1, MPI_COUNT_TYPE, commgrid->GetWorld());

    /*
     * Initialize displacement parameters now that we know both send and receive counts.
     */
    sdispls.front() = rdispls.front() = 0;

    std::partial_sum(sendcnt.begin(), sendcnt.end()-1, sdispls.begin()+1);
    std::partial_sum(recvcnt.begin(), recvcnt.end()-1, rdispls.begin()+1);

    /*
     * Need total send and receive buffer sizes in order to allocate
     * the memory.
     */
    int64_t totsend = std::accumulate(sendcnt.begin(), sendcnt.end(), 0);
    int64_t totrecv = std::accumulate(recvcnt.begin(), recvcnt.end(), 0);

    /*
     * Pack send buffer.
     */
    Vector<uint8_t> sendbuf(totsend, 0);

    for (int i = 0; i < nprocs; ++i)
    {
        assert(kmerbuckets[i].size() == (sendcnt[i] / TKmer::N_BYTES));

        /*
         * Get starting adddress of buffer space for kmerbuckets[i].
         */
        uint8_t *addrs2fill = sendbuf.data() + sdispls[i];

        for (MPI_Count_type j = 0; j < kmerbuckets[i].size(); ++j)
        {
            /*
             * Copy jth k-mer to be sent to ith processor into the
             * next available spot in the buffer.
             */
            (kmerbuckets[i][j]).CopyDataInto(addrs2fill);
            addrs2fill += TKmer::N_BYTES;
        }

        /*
         * Now that the k-mers to be sent to the ith processor
         * have been fully packed into the buffer, we don't
         * need them anymore.
         */
        kmerbuckets[i].clear();
    }

    /*
     * Allocate receive buffer.
     */
    Vector<uint8_t> recvbuf(totrecv, 0);

    /*
     * Send all the k-mers around.
     */
    MPI_ALLTOALLV(sendbuf.data(), sendcnt.data(), sdispls.data(), MPI_BYTE, recvbuf.data(), recvcnt.data(), rdispls.data(), MPI_BYTE, commgrid->GetWorld());

    /*
     * Get actual number of k-mer received.
     */
    uint64_t nkmersrecvd = totrecv / TKmer::N_BYTES;

    uint8_t *addrs2read = recvbuf.data();

    for (uint64_t i = 0; i < nkmersrecvd; ++i)
    {
        TKmer mer(addrs2read);
        addrs2read += TKmer::N_BYTES;

        if (kmermap.find(mer) == kmermap.end())
        {
            kmermap.insert({mer, KmerCountEntry()});
        }
    }

    return kmermap;
}

void GetKmerCountMapValues(const Vector<String>& myreads, KmerCountMap& kmermap, SharedPtr<CommGrid> commgrid)
{
    return;
}

