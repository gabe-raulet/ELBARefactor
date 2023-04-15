#ifndef KMER_COMM_H_
#define KMER_COMM_H_

#include "common.h"
#include "Kmer.h"
#include "CommGrid.h"

#ifndef LOWER_KMER_FREQ
#define LOWER_KMER_FREQ 20
#endif

#ifndef UPPER_KMER_FREQ
#define UPPER_KMER_FREQ 30
#endif

typedef uint16_t PosInRead;
typedef uint64_t ReadId;

typedef Array<PosInRead, UPPER_KMER_FREQ> POSITIONS;
typedef Array<ReadId,    UPPER_KMER_FREQ> READIDS;

typedef Tuple<READIDS, POSITIONS, int> KmerCountEntry;
typedef Map<TKmer, KmerCountEntry> KmerCountMap;

Set<TKmer> GetKmerCountMapKeys(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid);
void GetKmerCountMapValues(const Vector<String>& myreads, KmerCountMap& kmermap, SharedPtr<CommGrid> commgrid);

#endif
