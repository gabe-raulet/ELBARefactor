#ifndef TIMER_H_
#define TIMER_H_

#include <mpi.h>

class Timer
{
    MPI_Comm comm;
    double resettime;
    int myrank, nprocs;

public:
    Timer(MPI_Comm comm) : comm(comm)
    {
        MPI_Comm_rank(comm, &myrank);
        MPI_Comm_size(comm, &nprocs);
        Reset();
    }

    void Reset()
    {
        MPI_Barrier(comm);
        resettime = MPI_Wtime();
    }

    void GetTimes(double& mintime, double& maxtime, double& avgtime) const
    {
        double mytime = MPI_Wtime() - resettime;
        MPI_Allreduce(&mytime, &mintime, 1, MPI_DOUBLE, MPI_MIN, comm);
        MPI_Allreduce(&mytime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, comm);
        MPI_Allreduce(&mytime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, comm);
        avgtime /= nprocs;
    }
};

#endif
