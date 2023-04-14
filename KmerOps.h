#ifndef KMER_OPS_H_
#define KMER_OPS_H_

#include "common.h"
#include "Kmer.h"
#include "CommGrid.h"

#define MAX_ALLTOALL_MEM (128*1024*1024) /* 128 MB */

Set<TKmer> GetLocalKmers(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid);
Set<TKmer> GetKmersSmart(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid);

#endif
