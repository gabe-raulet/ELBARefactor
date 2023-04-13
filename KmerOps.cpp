#include "KmerOps.h"
#include "Bloom.h"
#include "Buffer.h"
#include <numeric>
#include <algorithm>

constexpr int64_t cardinality = 10000;

bool KmerParsePass
(
    Set<TKmer>&            localkmers, /* received k-mers will go here */
    const Vector<String>&  myreads,
    size_t&                myoffset,
    Bloom&                 bm,
    buffer_t              *scratch1,
    buffer_t              *scratch2,
    SharedPtr<CommGrid>    commgrid
)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    int len;
    size_t nreads = myreads.size();
    size_t maxsending = 0;
    size_t kmersthisbatch = 0;
    size_t bytesperkmer = TKmer::N_BYTES;
    size_t bytesperentry = bytesperkmer + 4;
    size_t memthreshold = (MAX_ALLTOALL_MEM / nprocs) * 2;

    Vector<TKmer> mykmers;
    Vector<Vector<TKmer>> outgoing(nprocs);

    while (myoffset < nreads)
    {
        String myseq = myreads[myoffset];
        len = myseq.size();

        if (len <= TKmer::k)
            continue;

        Vector<TKmer> kmers = TKmer::GetRepKmers(myseq);
        kmersthisbatch += kmers.size();

        size_t maxsending = 0;

        for (size_t j = 0; j < kmers.size(); ++j)
        {
            size_t sending;
            int owner = kmers[j].GetInferredOwner(nprocs);
            outgoing[owner].push_back(kmers[j]);
            sending = outgoing[owner].size();
            maxsending = std::max(sending, maxsending);
        }

        myoffset++;

        if (maxsending * bytesperentry >= memthreshold || (kmersthisbatch + len) * bytesperentry >= MAX_ALLTOALL_MEM)
            break;
    }

    bool finished = (myoffset >= nreads);

    std::vector<int> sendcnt(nprocs);
    std::vector<int> recvcnt(nprocs);
    std::vector<int> sdispls(nprocs);
    std::vector<int> rdispls(nprocs);

    for (int i = 0; i < nprocs; ++i)
    {
        sendcnt[i] = outgoing[i].size() * bytesperkmer;
    }

    MPI_Alltoall(sendcnt.data(), 1, MPI_INT, recvcnt.data(), 1, MPI_INT, commgrid->GetWorld());

    sdispls.front() = 0;
    rdispls.front() = 0;

    std::partial_sum(sendcnt.begin(), sendcnt.end()-1, sdispls.begin()+1);
    std::partial_sum(recvcnt.begin(), recvcnt.end()-1, rdispls.begin()+1);

    int64_t totsend = std::accumulate(sendcnt.begin(), sendcnt.end(), 0);
    int64_t totrecv = std::accumulate(recvcnt.begin(), recvcnt.end(), 0);

    grow_buffer(scratch1, totsend);

    uint8_t *sendbuf = (uint8_t*) get_start_buffer(scratch1);

    for (int i = 0; i < nprocs; ++i)
    {
        size_t nkmers2send = outgoing[i].size();
        uint8_t *addrs2fill = sendbuf + sdispls[i];

        for (size_t j = 0; j < nkmers2send; ++j)
        {
            (outgoing[i][j]).CopyDataInto(addrs2fill);
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

    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    Bloom bm(cardinality, 0.05);

    buffer_t *scratch1 = init_buffer(MAX_ALLTOALL_MEM);
    buffer_t *scratch2 = init_buffer(MAX_ALLTOALL_MEM);

    size_t myoffset = 0;
    size_t numreads = myreads.size();

    bool finished;

    do
    {
        finished = KmerParsePass(localkmers, myreads, myoffset, bm, scratch1, scratch2, commgrid);
    } while (!finished);

    free_buffer(scratch1);
    free_buffer(scratch2);

    return localkmers;
}
