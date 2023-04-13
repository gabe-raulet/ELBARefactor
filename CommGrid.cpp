#include "CommGrid.h"
#include <cmath>
#include <cassert>
#include <iostream>

CommGrid::CommGrid(MPI_Comm comm)
{
    int nprocs;

    MPI_Comm_dup(comm, &world);
    MPI_Comm_rank(world, &myrank);
    MPI_Comm_size(world, &nprocs);

    grcols = grrows = static_cast<int>(std::sqrt(static_cast<float>(nprocs)));

    if (grcols * grrows != nprocs)
    {
        std::cerr << "Must have a square logical processor grid" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    myproccol = myrank % grcols;
    myprocrow = myrank / grcols;

    /*
     * Create row and column communicators. MPI_Comm_split is called
     * collectively by all processes, and the processes called
     * with the same color are put in the same new communicator.
     */

    MPI_Comm_split(world, myprocrow, myrank, &row_world);
    MPI_Comm_split(world, myproccol, myrank, &col_world);

    /*
     * Sanity check. The rank of a process within its
     * row communicator is the same as the processes'
     * column within the grid, and vice versa for its
     * column communicator.
     */

    int row_rank, col_rank;
    MPI_Comm_rank(row_world, &row_rank);
    MPI_Comm_rank(col_world, &col_rank);

    assert(row_rank == myproccol);
    assert(col_rank == myprocrow);
}

CommGrid::CommGrid(const CommGrid& rhs) : grrows(rhs.grrows), grcols(rhs.grcols), myprocrow(rhs.myprocrow), myproccol(rhs.myproccol), myrank(rhs.myrank)
{
    MPI_Comm_dup(rhs.world, &world);
    MPI_Comm_dup(rhs.row_world, &row_world);
    MPI_Comm_dup(rhs.col_world, &col_world);
}

CommGrid::~CommGrid()
{
    MPI_Comm_free(&world);
    MPI_Comm_free(&row_world);
    MPI_Comm_free(&col_world);
}

CommGrid& CommGrid::operator=(const CommGrid& rhs)
{
    if (this != &rhs)
    {
        MPI_Comm_free(&world);
        MPI_Comm_free(&row_world);
        MPI_Comm_free(&col_world);

        grrows = rhs.grrows;
        grcols = rhs.grcols;
        myrank = rhs.myrank;
        myprocrow = rhs.myprocrow;
        myproccol = rhs.myproccol;

        MPI_Comm_dup(rhs.world, &world);
        MPI_Comm_dup(rhs.row_world, &row_world);
        MPI_Comm_dup(rhs.col_world, &col_world);
    }

    return *this;
}
