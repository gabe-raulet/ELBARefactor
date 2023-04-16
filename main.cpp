#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cassert>
#include <mpi.h>
#include "common.h"
#include "Kmer.h"
#include "KmerComm.h"
#include "CommGrid.h"
#include "FastaIndex.h"
#include "Timer.h"

Timer *timer = nullptr;
const String fasta_fname = "reads.fa";

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    {
        auto commgrid = SharedPtr<CommGrid>(new CommGrid(MPI_COMM_WORLD));
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();

        timer = new Timer(MPI_COMM_WORLD);

        FastaIndex index(fasta_fname, commgrid);

        Vector<String> myreads = FastaIndex::GetMyReads(index);

        KmerCountMap kmercounts = GetKmerCountMapKeys(myreads, commgrid);

        unsigned long numkmers = kmercounts.size();
        MPI_Allreduce(MPI_IN_PLACE, &numkmers, 1, MPI_UNSIGNED_LONG, MPI_SUM, commgrid->GetWorld());

        if (!myrank) std::cout << "A total of " << numkmers << " k-mers exist in the dataset" << std::endl;

        delete timer;
    }


    MPI_Finalize();
    return 0;
}

