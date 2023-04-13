#include "KmerOps.h"
#include "Bloom.h"

constexpr int64_t cardinality = 10000;

Set<TKmer> GetLocalKmers(const Vector<String>& myreads, SharedPtr<CommGrid> commgrid)
{
    Set<TKmer> localkmers;

    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    Bloom bm(cardinality, 0.05);



    return localkmers;
}
