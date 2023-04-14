#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cassert>
#include <mpi.h>
#include "common.h"
#include "Kmer.h"
#include "KmerOps.h"
#include "CommGrid.h"
#include "FastaIndex.h"

const int kmer_size = 31;
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

        Set<TKmer> localkmers = GetLocalKmers(myreads, commgrid);

        std::ostringstream ss;
        ss << "kmers.rank" << commgrid->GetRank() << ".txt";

        std::ofstream(ss.str());

        for (int i = 0; i < commgrid->GetSize(); ++i)
        {
            if (i == commgrid->GetRank())
            {
                for (auto itr = localkmers.begin(); itr != localkmers.end(); ++itr)
                    std::cout << *itr << std::endl;
            }
        }



        // Map<TKmer, KmerCountType> kmercounts;
        // KmerPass<false>(myreads, kmercounts, commgrid);
        // KmerPass<true>(myreads, kmercounts, commgrid);

        // uint64_t offset = myreads.size();
        // MPI_Exscan(MPI_IN_PLACE, &offset, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
        // if (commgrid->GetRank() == 0) offset = 0;

        // for (int i = 0; i < commgrid->GetSize(); ++i)
        // {
            // if (i == commgrid->GetRank())
            // {
                // for (auto itr = kmercounts.begin(); itr != kmercounts.end(); ++itr)
                // {
                    // for (int j = 0; j < UPPER_KMER_FREQ; ++j)
                    // {
                        // if (std::get<0>(itr->second)[j] != 0)
                            // std::cout << (itr->first).GetString() << "\t" << std::get<0>(itr->second)[j]-1+offset << "\t" << std::get<1>(itr->second)[j] << "\t" << std::get<2>(itr->second) << std::endl;
                    // }
                // }
            // }

            // MPI_Barrier(MPI_COMM_WORLD);

        // }

    }


    MPI_Finalize();
    return 0;
}

