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

const String fasta_fname = "reads.fa";

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    {
        auto commgrid = SharedPtr<CommGrid>(new CommGrid(MPI_COMM_WORLD));
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();

        if (!myrank)
        {
            std::cerr << "[ELBA]: -DKMER_SIZE=" << KMER_SIZE << " -DLOWER_KMER_FREQ=" << LOWER_KMER_FREQ << " -DUPPER_KMER_FREQ=" << UPPER_KMER_FREQ << " -DHLLFLAG=" << HLLFLAG << std::endl;
        }

        MPI_Barrier(commgrid->GetWorld());

        FastaIndex index(fasta_fname, commgrid);

        Vector<String> myreads = FastaIndex::GetMyReads(index);

        KmerCountMap kmercounts = GetKmerCountMapKeys(myreads, commgrid);

        unsigned long numkmers = kmercounts.size();
        MPI_Allreduce(MPI_IN_PLACE, &numkmers, 1, MPI_UNSIGNED_LONG, MPI_SUM, commgrid->GetWorld());

        if (!myrank) std::cout << "A total of " << numkmers << " k-mers exist in the dataset" << std::endl;

    }


    MPI_Finalize();
    return 0;
}

