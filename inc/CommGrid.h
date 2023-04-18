#ifndef COMM_GRID_H_
#define COMM_GRID_H_

#include <mpi.h>

class CommGrid
{
public:
    CommGrid(MPI_Comm comm);
    CommGrid(const CommGrid& rhs);
    ~CommGrid();

    CommGrid& operator=(const CommGrid& rhs);
    bool operator==(const CommGrid& rhs) const;
    bool operator!=(const CommGrid& rhs) const { return !(*this == rhs); }

    int GetRank(int rowrank, int colrank) { return rowrank * grcols + colrank; }
    int GetRank() { return myrank; }
    int GetRankInProcRow() { return myproccol; }
    int GetRankInProcCol() { return myprocrow; }
    int GetComplementRank() { return (grcols * myproccol) + myprocrow; }

    MPI_Comm& GetWorld() { return world; }
    MPI_Comm& GetRowWorld() { return row_world; }
    MPI_Comm& GetColWorld() { return col_world; }
    MPI_Comm GetWorld() const { return world; }
    MPI_Comm GetRowWorld() const { return row_world; }
    MPI_Comm GetColWorld() const { return col_world; }

    int GetGridRows() { return grrows; }
    int GetGridCols() { return grcols; }
    int GetSize() { return grrows * grcols; }

private:

    /*
     * MPI intracommunicators for ..
     */
    MPI_Comm world;     /* .. entire 2D processor grid */
    MPI_Comm row_world; /* .. grid rows */
    MPI_Comm col_world; /* .. grid columns */

    /*
     * Processor grid is (grrow X grcols)
     */
    int grrows, grcols;
    int myprocrow;
    int myproccol;
    int myrank;
};

#endif
