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

const int kmer_size = KMER_SIZE;
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

        KmerCountMap kmercounts = GetKmerCountMapKeys(myreads, commgrid);

        for (int i = 0; i < commgrid->GetSize(); ++i)
        {
            if (i == commgrid->GetRank())
            {
                for (auto itr = kmercounts.begin(); itr != kmercounts.end(); ++itr)
                {
                    std::cout << itr->first << std::endl;
                }
                std::cout << std::endl;
            }

            MPI_Barrier(commgrid->GetWorld());
        }

    }


    MPI_Finalize();
    return 0;
}

