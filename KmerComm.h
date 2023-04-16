#ifndef KMER_COMM_H_
#define KMER_COMM_H_

#include "common.h"
#include "Kmer.h"
#include "CommGrid.h"

#ifndef LOWER_KMER_FREQ
#error "LOWER_KMER_FREQ must be defined"
#endif

#ifndef UPPER_KMER_FREQ
#error "UPPER_KMER_FREQ must be defined"
#endif

#if (LOWER_KMER_FREQ > UPPER_KMER_FREQ)
#error "LOWER_KMER_FREQ must be less than or equal to UPPER_KMER_FREQ"
#endif

typedef uint16_t PosInRead;
typedef uint64_t ReadId;

typedef Array<PosInRead, UPPER_KMER_FREQ> POSITIONS;
typedef Array<ReadId,    UPPER_KMER_FREQ> READIDS;

typedef Tuple<READIDS, POSITIONS, int> KmerCountEntry;
typedef Map<TKmer, KmerCountEntry> KmerCountMap;

KmerCountMap GetKmerCountMapKeys(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid);
void GetKmerCountMapValues(const Vector<String>& myreads, KmerCountMap& kmermap, SharedPtr<CommGrid> commgrid);

struct KmerEstimateHandler
{
    HyperLogLog& hll;

    KmerEstimateHandler(HyperLogLog& hll) : hll(hll) {}

    void operator()(const TKmer& kmer)
    {
        hll.Add(kmer.GetString().c_str());
    }
};

struct KmerPartitionHandler
{
    int nprocs;
    Vector<Vector<TKmer>>& kmerbuckets;

    KmerPartitionHandler(Vector<Vector<TKmer>>& kmerbuckets) : nprocs(kmerbuckets.size()), kmerbuckets(kmerbuckets) {}

    void operator()(const TKmer& kmer)
    {
        uint64_t myhash = kmer.GetHash();
        double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
        size_t owner = range / std::numeric_limits<uint64_t>::max();
        assert(owner >= 0 && owner < static_cast<int>(nprocs));
        kmerbuckets[owner].push_back(kmer);
    }
};

template <typename KmerHandler>
void ForeachKmer(const Vector<String>& myreads, KmerHandler& handler)
{
    /*
     * Go through each local read.
     */
    for (auto readitr = myreads.begin(); readitr != myreads.end(); ++readitr)
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

        /*
         * Go through each k-mer seed.
         */
        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr)
        {
            handler(*meritr);
        }
    }
}

#endif
