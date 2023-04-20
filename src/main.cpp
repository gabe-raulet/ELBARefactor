#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cassert>
#include <numeric>
#include <mpi.h>
#include "common.h"
#include "Kmer.h"
#include "KmerComm.h"
#include "FastaIndex.h"
#include "ReadOverlap.h"
#include "KmerIntersect.h"

String fasta_fname = "data/reads.fa";

void LogAll(const String mylog, SharedPtr<CommGrid> commgrid);
String ProcessorName(SharedPtr<CommGrid> commgrid);
void PrintKmerHistogram(const KmerCountMap& kmermap, SharedPtr<CommGrid>& commgrid);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if (argc == 2) fasta_fname.assign(argv[1]);

    {
        MPI_Comm gridworld = MPI_COMM_WORLD;
        auto commgrid = SharedPtr<CommGrid>(new CommGrid(gridworld, 0, 0));
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();

        TuplePair<std::ostringstream> logstream;
        FastaIndex index(fasta_fname, commgrid);
        Vector<String> myreads = FastaIndex::GetMyReads(index);
        KmerCountMap kmercounts = GetKmerCountMapKeys(myreads, commgrid);

        size_t numkmers = kmercounts.size();
        MPI_Allreduce(MPI_IN_PLACE, &numkmers, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());

        if (!myrank) std::cout << "A total of " << numkmers << " k-mers exist in the dataset" << std::endl;

        GetKmerCountMapValues(myreads, kmercounts, commgrid);

        auto itr = kmercounts.begin();

        while (itr != kmercounts.end())
        {
            itr = std::get<2>(itr->second) < LOWER_KMER_FREQ? kmercounts.erase(itr) : ++itr;
        }

        PrintKmerHistogram(kmercounts, commgrid);

        uint64_t kmerid = kmercounts.size();
        uint64_t totkmers = kmerid;
        uint64_t totreads = myreads.size();

        MPI_Allreduce(&kmerid,      &totkmers, 1, MPI_UINT64_T, MPI_SUM, commgrid->GetWorld());
        MPI_Allreduce(MPI_IN_PLACE, &totreads, 1, MPI_UINT64_T, MPI_SUM, commgrid->GetWorld());

        MPI_Exscan(MPI_IN_PLACE, &kmerid, 1, MPI_UINT64_T, MPI_SUM, commgrid->GetWorld());
        if (myrank == 0) kmerid = 0;

        Vector<uint64_t> local_rowids, local_colids;
        Vector<PosInRead> local_positions;

        for (auto itr = kmercounts.begin(); itr != kmercounts.end(); ++itr)
        {
            READIDS& readids = std::get<0>(itr->second);
            POSITIONS& positions = std::get<1>(itr->second);
            int cnt = std::get<2>(itr->second);

            for (int j = 0; j < cnt; ++j)
            {
                local_colids.push_back(kmerid);
                local_rowids.push_back(readids[j]);
                local_positions.push_back(positions[j]);
            }

            kmerid++;
        }

        CT<uint64_t>::PDistVec drows(local_rowids, commgrid);
        CT<uint64_t>::PDistVec dcols(local_colids, commgrid);
        CT<PosInRead>::PDistVec dvals(local_positions, commgrid);

        CT<PosInRead>::PSpParMat A(totreads, totkmers, drows, dcols, dvals, true);

        auto AT = A;
        AT.Transpose();

        CT<ReadOverlap>::PSpParMat B = Mult_AnXBn_DoubleBuff<KmerIntersect, ReadOverlap, CT<ReadOverlap>::PSpDCCols>(A, AT);

        B.Prune([](const auto& item) { return item.count <= 1; });

        B.ParallelWriteMM("B.mtx", false, OverlapHandler());
    }

    MPI_Finalize();
    return 0;
}

void PrintKmerHistogram(const KmerCountMap& kmermap, SharedPtr<CommGrid>& commgrid)
{
    int maxcount = std::accumulate(kmermap.cbegin(), kmermap.cend(), 0, [](int cur, const auto& entry) { return std::max(cur, std::get<2>(entry.second)); });

    MPI_Allreduce(MPI_IN_PLACE, &maxcount, 1, MPI_INT, MPI_MAX, commgrid->GetWorld());

    Vector<int> histo(maxcount+1, 0);

    for (auto itr = kmermap.cbegin(); itr != kmermap.cend(); ++itr)
    {
        int cnt = std::get<2>(itr->second);
        assert(cnt >= 1);
        histo[cnt]++;
    }

    MPI_Allreduce(MPI_IN_PLACE, histo.data(), maxcount+1, MPI_INT, MPI_SUM, commgrid->GetWorld());

    int myrank = commgrid->GetRank();

    if (!myrank)
    {
        std::cout << "#count\tnumkmers" << std::endl;

        for (int i = 1; i < histo.size(); ++i)
        {
            if (histo[i] > 0)
            {
                std::cout << i << "\t" << histo[i] << std::endl;
            }
        }
    }
}

void LogAll(const String mylog, SharedPtr<CommGrid> commgrid)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    Vector<int> recvcnt, displs;
    int sendcnt = mylog.size();
    int totrecv;

    if (!myrank) recvcnt.resize(nprocs);

    MPI_Gather(&sendcnt, 1, MPI_INT, recvcnt.data(), 1, MPI_INT, 0, comm);

    if (!myrank)
    {
        displs.resize(nprocs);
        displs.front() = 0;
        std::partial_sum(recvcnt.begin(), recvcnt.end()-1, displs.begin()+1);
        totrecv = recvcnt.back() + displs.back();
    }

    Vector<char> recvbuf(totrecv);

    MPI_Gatherv(mylog.c_str(), sendcnt, MPI_CHAR, recvbuf.data(), recvcnt.data(), displs.data(), MPI_CHAR, 0, comm);

    if (!myrank)
    {
        char const *buf = recvbuf.data();

        for (int i = 0; i < nprocs; ++i)
        {
            String message(buf + displs[i], recvcnt[i]);
            std::cout << "processor[" << i+1 << "/" << nprocs << "] says '" << message << "'" << std::endl;
        }
    }

    MPI_Barrier(comm);
}

String ProcessorName(SharedPtr<CommGrid> commgrid)
{
    static bool initialized = false;
    static String name;

    if (!initialized)
    {
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();
        int rowrank = commgrid->GetRankInProcCol();
        int colrank = commgrid->GetRankInProcRow();
        std::ostringstream ss;
        ss  << "processor[" << myrank+1 << "/" << nprocs << "]...grid[" << rowrank+1 <<"," << colrank+1 << "])";
        name = ss.str();
        initialized = true;
    }

    return name;
}
