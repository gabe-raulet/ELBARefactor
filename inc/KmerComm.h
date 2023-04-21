#ifndef KMER_COMM_H_
#define KMER_COMM_H_

#include "common.h"
#include "Kmer.h"

#ifndef LOWER_KMER_FREQ
#error "LOWER_KMER_FREQ must be defined"
#endif

#ifndef UPPER_KMER_FREQ
#error "UPPER_KMER_FREQ must be defined"
#endif

#if (LOWER_KMER_FREQ > UPPER_KMER_FREQ)
#error "LOWER_KMER_FREQ must be less than or equal to UPPER_KMER_FREQ"
#endif

#ifndef MAX_ALLTOALL_MEM
// #define MAX_ALLTOALL_MEM (1024)
// #define MAX_ALLTOALL_MEM (128 * 1024)
#define MAX_ALLTOALL_MEM (128 * 1024 * 1024)
// #define MAX_ALLTOALL_MEM (32 * 1024 * 1024)
#endif

typedef uint16_t PosInRead;
typedef uint64_t ReadId;

typedef Array<PosInRead, UPPER_KMER_FREQ> POSITIONS;
typedef Array<ReadId,    UPPER_KMER_FREQ> READIDS;

typedef Tuple<TKmer, ReadId, PosInRead> KmerSeed;
typedef Tuple<READIDS, POSITIONS, int> KmerCountEntry;
typedef Map<TKmer, KmerCountEntry> KmerCountMap;

KmerCountMap GetKmerCountMapKeys(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid);
void GetKmerCountMapValues(const Vector<String>& myreads, KmerCountMap& kmermap, SharedPtr<CommGrid> commgrid);
int GetKmerOwner(const TKmer& kmer, int nprocs);

struct BatchState
{
    SharedPtr<CommGrid> commgrid;
    size_t mynumreads;
    size_t memthreshold;
    size_t mykmerssofar;
    size_t mymaxsending;
    ReadId myreadid;

    BatchState(size_t mynumreads, SharedPtr<CommGrid> commgrid) : commgrid(commgrid), mynumreads(mynumreads), memthreshold((MAX_ALLTOALL_MEM / commgrid->GetSize()) << 1), mykmerssofar(0), mymaxsending(0), myreadid(0) {}

    bool ReachedThreshold(const size_t len)
    {
        return (mymaxsending * TKmer::N_BYTES >= memthreshold || (mykmerssofar + len) * TKmer::N_BYTES >= MAX_ALLTOALL_MEM);
    }

    bool Finished() const
    {
        int alldone;
        int imdone = !!(myreadid >= static_cast<ReadId>(mynumreads));
        MPI_Allreduce(&imdone, &alldone, 1, MPI_INT, MPI_LAND, commgrid->GetWorld());
        return bool(alldone);
    }
};

struct KmerEstimateHandler
{
    HyperLogLog& hll;

    KmerEstimateHandler(HyperLogLog& hll) : hll(hll) {}

    void operator()(const TKmer& kmer, size_t kid, size_t rid)
    {
        hll.Add(kmer.GetString().c_str());
    }
};

struct KmerPartitionHandler
{
    int nprocs;
    Vector<Vector<TKmer>>& kmerbuckets;

    KmerPartitionHandler(Vector<Vector<TKmer>>& kmerbuckets) : nprocs(kmerbuckets.size()), kmerbuckets(kmerbuckets) {}

    void operator()(const TKmer& kmer, size_t kid, size_t rid)
    {
        kmerbuckets[GetKmerOwner(kmer, nprocs)].push_back(kmer);
    }

    void operator()(const TKmer& kmer, BatchState& state, size_t kid)
    {
        auto& kmerbucket = kmerbuckets[GetKmerOwner(kmer, nprocs)];
        kmerbucket.push_back(kmer);
        state.mymaxsending = std::max(kmerbucket.size(), state.mymaxsending);
    }
};

struct KmerParserHandler
{
    int nprocs;
    ReadId readoffset;
    Vector<Vector<KmerSeed>>& kmerseeds;

    KmerParserHandler(Vector<Vector<KmerSeed>>& kmerseeds, ReadId readoffset) : nprocs(kmerseeds.size()), readoffset(readoffset), kmerseeds(kmerseeds) {}

    void operator()(const TKmer& kmer, size_t kid, size_t rid)
    {
        kmerseeds[GetKmerOwner(kmer, nprocs)].emplace_back(kmer, static_cast<ReadId>(rid) + readoffset, static_cast<PosInRead>(kid));
    }

};

template <typename KmerHandler>
void ForeachKmer(const Vector<String>& myreads, KmerHandler& handler)
{
    size_t i = 0;
    /*
     * Go through each local read.
     */
    for (auto readitr = myreads.begin(); readitr != myreads.end(); ++readitr, ++i)
    {
        /*
         * If it is too small then continue to the next one.
         */
        if (readitr->size() < KMER_SIZE)
            continue;

        /*
         * Get all the representative k-mer seeds.
         */
        Vector<TKmer> repmers = TKmer::GetRepKmers(*readitr);

        size_t j = 0;

        /*
         * Go through each k-mer seed.
         */
        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr, ++j)
        {
            handler(*meritr, j, i);
        }
    }
}


template <typename KmerHandler>
void ForeachKmer(const Vector<String>& myreads, KmerHandler& handler, BatchState& state)
{
    for (; state.myreadid < static_cast<ReadId>(myreads.size()); state.myreadid++)
    {
        const String& sequence = myreads[state.myreadid];

        if (sequence.size() < KMER_SIZE)
            continue;

        Vector<TKmer> repmers = TKmer::GetRepKmers(sequence);
        state.mykmerssofar += repmers.size();

        size_t j = 0;

        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr, ++j)
        {
            handler(*meritr, state, j);
        }

        if (state.ReachedThreshold(sequence.size()))
        {
            state.myreadid++;
            return;
        }
    }
}

#endif
