#ifndef KMER_OPS_H_
#define KMER_OPS_H_

#include "common.h"
#include "Kmer.h"
#include "CommGrid.h"
#include <cstdint>

#ifndef LOWER_KMER_FREQ
#define LOWER_KMER_FREQ 20
#endif

#ifndef UPPER_KMER_FREQ
#define UPPER_KMER_FREQ 30
#endif

#ifndef MAX_ALLTOALL_MEM
#define MAX_ALLTOALL_MEM (128*1024*1024) /* 128 MB */
#endif

typedef uint16_t PosInRead;
typedef uint64_t ReadId;

typedef Array<PosInRead, UPPER_KMER_FREQ> POSITIONS;
typedef Array<ReadId,    UPPER_KMER_FREQ> READIDS;

typedef Tuple<READIDS, POSITIONS, int> KmerCountType;
typedef Tuple<TKmer, KmerCountType> KmerValue;

Set<TKmer> GetLocalKmers(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid);
Set<TKmer> GetKmersSmart(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid);

template <bool COUNT_PHASE>
void KmerPass(const Vector<String>& myreads, Map<TKmer, KmerCountType>& kmercounts, SharedPtr<CommGrid> commgrid);

#endif
