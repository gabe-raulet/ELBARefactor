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

        std::ostringstream ss;
        ss << "kmers.rank" << commgrid->GetRank() << ".txt";

        std::ofstream filestream(ss.str());

        for (auto itr = mykmers.begin(); itr != mykmers.end(); ++itr)
        {
            filestream << itr->GetString() << std::endl;
        }

        filestream.close();
    }


    MPI_Finalize();
    return 0;
}

