#include <iostream>
#include <sstream>
#include <cstdint>
#include <cassert>
#include <mpi.h>
#include "common.h"
#include "Kmer.h"
#include "CommGrid.h"
#include "FastaIndex.h"

const int kmer_size = 15;
const String fasta_fname = "reads.fa";

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    TKmer::SetKmerSize(kmer_size);

    {
        auto commgrid = SharedPtr<CommGrid>(new CommGrid(MPI_COMM_WORLD));

        FastaIndex index(fasta_fname, commgrid);

        index.PrintInfo();

        Vector<String> myreads = FastaIndex::GetMyReads(index);

        // Set<TKmer> mykmers = GetLocalKmers(myreads, commgrid);
    }


    MPI_Finalize();
    return 0;
}

