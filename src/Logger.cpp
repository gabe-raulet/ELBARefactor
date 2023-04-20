#include "Logger.h"

void LogAll(const String mylog, SharedPtr<CommGrid> commgrid)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    Vector<int> recvcnt, displs;
    Vector<char> recvbuf;
    int sendcnt = mylog.size();
    int totrecv;

    if (!myrank) recvcnt.resize(nprocs);

    MPI_Gather(&sendcnt, 1, MPI_INT, recvcnt.data(), 1, MPI_INT, 0, comm);

    if (!myrank)
    {
        displs.resize(nprocs);
        displs.front() = 0;
        std::partial_sum(recvcnt.begin(), recvcnt.end()-1, displs.begin()+1);
        recvbuf.resize(recvcnt.back() + displs.back());
    }

    MPI_Gatherv(mylog.c_str(), sendcnt, MPI_CHAR, recvbuf.data(), recvcnt.data(), displs.data(), MPI_CHAR, 0, comm);

    if (!myrank)
    {
        char const *buf = recvbuf.data();

        for (int i = 0; i < nprocs; ++i)
        {
            String message(buf + displs[i], recvcnt[i]);
            std::cout << "processor[" << i+1 << "/" << nprocs << "] says '" << message << "'" << std::endl;
        }
        std::cout << std::endl;
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
