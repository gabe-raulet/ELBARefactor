#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cassert>
#include <numeric>
#include <mpi.h>
#include "common.h"
#include "FastaIndex.h"
#include "Logger.h"

String fasta_fname = "data/reads.fa";

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if (argc == 2) fasta_fname.assign(argv[1]);

    {
        MPI_Comm gridworld = MPI_COMM_WORLD;
        auto commgrid = SharedPtr<CommGrid>(new CommGrid(gridworld, 0, 0));
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();

        if (!myrank)
        {
            std::cout << "-DKMER_SIZE=" << KMER_SIZE << " "
                      << "-DLOWER_KMER_FREQ=" << LOWER_KMER_FREQ << " "
                      << "-DUPPER_KMER_FREQ=" << UPPER_KMER_FREQ << " "
                      << "-DUSE_BLOOM=" << USE_BLOOM << "\n" << std::endl;
        }

        MPI_Barrier(gridworld);

        FastaIndex index(fasta_fname, commgrid);
        Vector<String> myreads = index.GetMyReads();

    }

    MPI_Finalize();
    return 0;
}
