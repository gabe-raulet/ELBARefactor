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

struct KmerEstimateHandler
{
    HyperLogLog& hll;

    KmerEstimateHandler(HyperLogLog& hll) : hll(hll) {}

    void operator()(const TKmer& kmer, const String& myread, size_t kid, size_t rid)
    {
        hll.Add(kmer.GetString().c_str());
    }
};

struct KmerPartitionHandler
{
    int nprocs;
    Vector<Vector<TKmer>>& kmerbuckets;

    KmerPartitionHandler(Vector<Vector<TKmer>>& kmerbuckets) : nprocs(kmerbuckets.size()), kmerbuckets(kmerbuckets) {}

    void operator()(const TKmer& kmer, const String& myread, size_t kid, size_t rid)
    {
        kmerbuckets[GetKmerOwner(kmer, nprocs)].push_back(kmer);
    }
};

struct KmerParserHandler
{
    int nprocs;
    ReadId readoffset;
    Vector<Vector<KmerSeed>>& kmerseeds;

    KmerParserHandler(Vector<Vector<KmerSeed>>& kmerseeds, ReadId readoffset) : nprocs(kmerseeds.size()), readoffset(readoffset), kmerseeds(kmerseeds) {}

    void operator()(const TKmer& kmer, const String& myread, size_t kid, size_t rid)
    {
        kmerseeds[GetKmerOwner(kmer, nprocs)].emplace_back(kmer, static_cast<ReadId>(rid + readoffset), static_cast<PosInRead>(kid));
    }

};

template <typename KmerHandler>
void ForeachKmer(const Vector<String>& myreads, KmerHandler& handler)
{
    size_t i = 0, j = 0;
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

        /*
         * Go through each k-mer seed.
         */
        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr, ++j)
        {
            handler(*meritr, *readitr, j, i);
        }
    }
}

#endif
