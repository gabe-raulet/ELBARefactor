#include <iostream>
#include <sstream>
#include <cstdint>
#include <cassert>
#include <mpi.h>
#include "common.h"
#include "Kmer.h"
#include "KmerOps.h"
#include "CommGrid.h"
#include "FastaIndex.h"

const int kmer_size = 7;
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

        Set<TKmer> mykmers = GetLocalKmers(myreads, commgrid);


        for (int i = 0; i < commgrid->GetSize(); ++i)
        {
            if (i == commgrid->GetRank())
                std::cerr << "Process " << i << " will store " << mykmers.size() << " k-mers" << std::endl;

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }


    MPI_Finalize();
    return 0;
}

