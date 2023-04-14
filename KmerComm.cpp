#include "KmerComm.h"

KmerCountMap GetKmerCountMapKeys(const Vector <String>& myreads, SharedPtr<CommGrid> commgrid)
{
    KmerCountMap kmermap;

    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    size_t numreads = myreads.size();


    return kmermap;
}

void GetKmerCountMapValues(const Vector<String>& myreads, KmerCountMap& kmermap, SharedPtr<CommGrid> commgrid)
{
    return;
}

