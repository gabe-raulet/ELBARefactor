#ifndef MPI_PRINT_H_
#define MPI_PRINT_H_

#include "common.h"
#include "CommGrid.h"
#include <cassert>
#include <limits>

template <class T>
void PrintVector(const Vector<T>& vec, SharedPtr<CommGrid>& commgrid)
{
    MPI_Comm world = commgrid->GetWorld();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    assert(vec.size() <= Limits<int>::max());

    for (int i = 0; i < nprocs; ++i)
    {
        if (myrank == i)
        {
            std::cout << myrank << ": ";
            int len = vec.size();

            for (int j = 0; j < len-1; ++j)
            {
                std::cout << vec[j] << ", ";
            }

            std::cout << vec[len-1] << std::endl;
        }

        MPI_Barrier(world);
    }
}

template <class T>
void PrintItem(const T& item, SharedPtr<CommGrid>& commgrid)
{
    Vector<T> v = {item};
    PrintVector<T>(v, commgrid);
}

template <class T>
void PrintItem(const Vector<T>& vec, SharedPtr<CommGrid>& commgrid)
{
    PrintVector<T>(vec, commgrid);
}



#endif
