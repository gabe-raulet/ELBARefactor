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

String fasta_fname = "data/reads.fa";

template <bool HAS_SEEDS>
void PrintKmerMap(KmerCountMap& kmermap, SharedPtr<CommGrid>& commgrid);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if (argc == 2) fasta_fname.assign(argv[1]);

    {
        auto commgrid = SharedPtr<CommGrid>(new CommGrid(MPI_COMM_WORLD));
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();

        FastaIndex index(fasta_fname, commgrid);

        Vector<String> myreads = FastaIndex::GetMyReads(index);

        KmerCountMap kmercounts = GetKmerCountMapKeys(myreads, commgrid);

        PrintKmerMap<false>(kmercounts, commgrid);

        size_t numkmers = kmercounts.size();
        MPI_Allreduce(MPI_IN_PLACE, &numkmers, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());

        if (!myrank) std::cout << "A total of " << numkmers << " k-mers exist in the dataset" << std::endl;

        GetKmerCountMapValues(myreads, kmercounts, commgrid);

        PrintKmerMap<true>(kmercounts, commgrid);

    }


    MPI_Finalize();
    return 0;
}

template <bool HAS_SEEDS>
void PrintKmerMap(KmerCountMap& kmermap, SharedPtr<CommGrid>& commgrid)
{
    int myrank = commgrid->GetRank();

    std::ostringstream ss;

    if (HAS_SEEDS)
    {
        ss << "output.seeds.rank" << myrank << ".txt";
    }
    else
    {
        ss << "output.kmers.rank" << myrank << ".txt";
    }

    std::ofstream filestream(ss.str());

    for (auto itr = kmermap.begin(); itr != kmermap.end(); ++itr)
    {
        String kmer = (itr->first).GetString();

        if (HAS_SEEDS)
        {
            READIDS& readids = std::get<0>(itr->second);
            POSITIONS& positions = std::get<1>(itr->second);
            int& count = std::get<2>(itr->second);

            for (int i = 0; i < count; ++i)
            {
                filestream << kmer << "\t" << readids[i] << "\t" << positions[i] << std::endl;
            }
        }
        else
        {
            filestream << kmer << std::endl;
        }
    }

    filestream.close();
}
