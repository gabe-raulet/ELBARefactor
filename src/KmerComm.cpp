#include "KmerComm.h"
#include "Bloom.h"
#include "Logger.h"
#include <cstring>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <cmath>

#if USE_BLOOM == 1
static Bloom *bm = nullptr;
#else
static_assert(USE_BLOOM == 0);
#endif

KmerCountMap GetKmerCountMapKeys(const Vector <String>& myreads, SharedPtr<CommGrid> commgrid)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    KmerCountMap kmermap;
    HyperLogLog hll;
    size_t numreads;
    size_t avgcardinality;
    double cardinality;
    double mycardinality;
    Vector<Vector<TKmer>> kmerbuckets(nprocs);
    Vector<MPI_Count_type> sendcnt(nprocs), recvcnt(nprocs);
    Vector<MPI_Displ_type> sdispls(nprocs), rdispls(nprocs);
    Vector<uint8_t> sendbuf, recvbuf;
    size_t totsend, totrecv;
    size_t numkmerseeds;

    numreads = myreads.size();
    KmerEstimateHandler estimator(hll);
    ForeachKmer(myreads, estimator);
    mycardinality = hll.Estimate();
    hll.ParallelMerge(commgrid->GetWorld());
    cardinality = hll.Estimate();
    avgcardinality = static_cast<size_t>(std::ceil(cardinality / nprocs));
    kmermap.reserve(avgcardinality);
    bm = new Bloom(static_cast<int64_t>(std::ceil(cardinality)), 0.05);
    KmerPartitionHandler partitioner(kmerbuckets);
    ForeachKmer(myreads, partitioner);
    for (int i = 0; i < nprocs; ++i)
    {
        sendcnt[i] = kmerbuckets[i].size() * TKmer::N_BYTES;
    }
    MPI_ALLTOALL(sendcnt.data(), 1, MPI_COUNT_TYPE, recvcnt.data(), 1, MPI_COUNT_TYPE, commgrid->GetWorld());
    sdispls.front() = rdispls.front() = 0;
    std::partial_sum(sendcnt.begin(), sendcnt.end()-1, sdispls.begin()+1);
    std::partial_sum(recvcnt.begin(), recvcnt.end()-1, rdispls.begin()+1);
    totsend = std::accumulate(sendcnt.begin(), sendcnt.end(), 0);
    totrecv = std::accumulate(recvcnt.begin(), recvcnt.end(), 0);
    sendbuf.resize(totsend, 0);
    for (int i = 0; i < nprocs; ++i)
    {
        uint8_t *dest = sendbuf.data() + sdispls[i];
        for (MPI_Count_type j = 0; j < kmerbuckets[i].size(); ++j)
        {
            memcpy(dest, kmerbuckets[i][j].GetBytes(), TKmer::N_BYTES);
            dest += TKmer::N_BYTES;
        }
        kmerbuckets[i].clear();
    }
    recvbuf.resize(totrecv, 0);
    MPI_ALLTOALLV(sendbuf.data(), sendcnt.data(), sdispls.data(), MPI_BYTE, recvbuf.data(), recvcnt.data(), rdispls.data(), MPI_BYTE, commgrid->GetWorld());
    numkmerseeds = totrecv / TKmer::N_BYTES;
    uint8_t *src = recvbuf.data();
    for (size_t i = 0; i < numkmerseeds; ++i)
    {
        TKmer mer(src);
        src += TKmer::N_BYTES;
        if (bm->Check(mer.GetBytes(), TKmer::N_BYTES))
        {
            if (kmermap.find(mer) == kmermap.end())
                kmermap.insert({mer, KmerCountEntry({}, {}, 0)});
        }
        else bm->Add(mer.GetBytes(), TKmer::N_BYTES);
    }
    return kmermap;
}

void GetKmerCountMapValues(const Vector<String>& myreads, KmerCountMap& kmermap, SharedPtr<CommGrid> commgrid)
{
    std::unique_ptr<std::ostringstream> logstream;
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    size_t numreads = myreads.size();
    Vector<Vector<KmerSeed>> kmerseeds(nprocs);
    size_t readoffset = numreads;
    MPI_Exscan(&numreads, &readoffset, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());
    if (!myrank) readoffset = 0;
    KmerParserHandler parser(kmerseeds, static_cast<ReadId>(readoffset));
    ForeachKmer(myreads, parser);
    Vector<MPI_Count_type> sendcnt(nprocs), recvcnt(nprocs);
    Vector<MPI_Displ_type> sdispls(nprocs), rdispls(nprocs);
    constexpr size_t seedbytes = TKmer::N_BYTES + sizeof(ReadId) + sizeof(PosInRead);
    logstream.reset(new std::ostringstream());
    *logstream << std::setprecision(4) << "sending 'row' k-mers to each processor in this amount (megabytes): {";
    for (int i = 0; i < nprocs; ++i)
    {
        sendcnt[i] = kmerseeds[i].size() * seedbytes;
        *logstream << (static_cast<double>(sendcnt[i]) / (1024 * 1024)) << ",";
    }
    *logstream << "}";
    LogAll(logstream->str(), commgrid);
    MPI_ALLTOALL(sendcnt.data(), 1, MPI_COUNT_TYPE, recvcnt.data(), 1, MPI_COUNT_TYPE, commgrid->GetWorld());
    sdispls.front() = rdispls.front() = 0;
    std::partial_sum(sendcnt.begin(), sendcnt.end()-1, sdispls.begin()+1);
    std::partial_sum(recvcnt.begin(), recvcnt.end()-1, rdispls.begin()+1);
    size_t totsend = std::accumulate(sendcnt.begin(), sendcnt.end(), 0);
    size_t totrecv = std::accumulate(recvcnt.begin(), recvcnt.end(), 0);
    Vector<uint8_t> sendbuf(totsend, 0);
    for (int i = 0; i < nprocs; ++i)
    {
        assert(kmerseeds[i].size() == (sendcnt[i] / seedbytes));
        uint8_t *addrs2fill = sendbuf.data() + sdispls[i];
        for (MPI_Count_type j = 0; j < kmerseeds[i].size(); ++j)
        {
            auto& seeditr = kmerseeds[i][j];
            TKmer kmer = std::get<0>(seeditr);
            ReadId readid = std::get<1>(seeditr);
            PosInRead pos = std::get<2>(seeditr);
            memcpy(addrs2fill, kmer.GetBytes(), TKmer::N_BYTES);
            memcpy(addrs2fill + TKmer::N_BYTES, &readid, sizeof(ReadId));
            memcpy(addrs2fill + TKmer::N_BYTES + sizeof(ReadId), &pos, sizeof(PosInRead));
            addrs2fill += seedbytes;
        }
        kmerseeds[i].clear();
    }
    Vector<uint8_t> recvbuf(totrecv, 0);
    MPI_ALLTOALLV(sendbuf.data(), sendcnt.data(), sdispls.data(), MPI_BYTE, recvbuf.data(), recvcnt.data(), rdispls.data(), MPI_BYTE, commgrid->GetWorld());
    size_t numkmerseeds = totrecv / seedbytes;
    logstream.reset(new std::ostringstream());
    *logstream << "received a total of " << numkmerseeds << " 'row' k-mers in second ALLTOALL exchange";
    LogAll(logstream->str(), commgrid);
    uint8_t *addrs2read = recvbuf.data();
    for (size_t i = 0; i < numkmerseeds; ++i)
    {
        TKmer kmer(addrs2read);
#if USE_BLOOM == 1
        if (!bm->Check(kmer.GetBytes(), TKmer::N_BYTES))
            continue;
#else
        static_assert(USE_BLOOM == 0);
#endif
        ReadId readid = *((ReadId*)(addrs2read + TKmer::N_BYTES));
        PosInRead pos = *((PosInRead*)(addrs2read + TKmer::N_BYTES + sizeof(ReadId)));
        addrs2read += seedbytes;
        auto kmitr = kmermap.find(kmer);
        if (kmitr == kmermap.end())
            continue;
        KmerCountEntry& entry = kmitr->second;
        READIDS& readids      = std::get<0>(entry);
        POSITIONS& positions  = std::get<1>(entry);
        int& count            = std::get<2>(entry);
        if (count >= UPPER_KMER_FREQ) /* TODO: There is probably a more efficient solution: deleting k-mer from kmermap (?) */
        {
            kmermap.erase(kmer);
            continue;
        }
        readids[count] = readid;
        positions[count] = pos;
        count++;
    }
    logstream.reset(new std::ostringstream());
    *logstream << numkmerseeds;
#if USE_BLOOM == 1
    *logstream << " row k-mers filtered by Bloom filter, hash table, and upper k-mer bound threshold into " << kmermap.size() << " semi-reliable 'column' k-mers";
    delete bm;
#else
    static_assert(USE_BLOOM == 0);
    *logstream << " row k-mers filtered by hash table and upper k-mer bound threshold into " << kmermap.size() << " semi-reliable 'column' k-mers";
#endif
    LogAll(logstream->str(), commgrid);
    auto itr = kmermap.begin();
    while (itr != kmermap.end())
    {
        itr = std::get<2>(itr->second) < LOWER_KMER_FREQ? kmermap.erase(itr) : ++itr;
    }
    size_t numkmers = kmermap.size();
    MPI_Allreduce(MPI_IN_PLACE, &numkmers, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());
    if (!myrank)
    {
        std::cout << "A total of " << numkmers << " reliable 'column' k-mers found\n" << std::endl;
    }
    MPI_Barrier(commgrid->GetWorld());
}

int GetKmerOwner(const TKmer& kmer, int nprocs)
{
    uint64_t myhash = kmer.GetHash();
    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / std::numeric_limits<uint64_t>::max();
    assert(owner >= 0 && owner < static_cast<int>(nprocs));
    return static_cast<int>(owner);
}
