#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cassert>
#include <cstdio>
#include <numeric>
#include <unistd.h>
#include <mpi.h>
#include "common.h"
#include "compiletime.h"
#include "FastaIndex.h"
#include "Logger.h"

int returncode;
String fasta_fname;

/*
 * MPI communicator info.
 */
MPI_Comm comm;
int myrank;
int nprocs;

/*
 * X-Drop alignment parameters.
 */
int mat = 1;           /* match score */
int mis = -1;          /* mismatch score */
int gap = -1;          /* gap score */
int xdrop_cutoff = 15; /* x-drop cutoff score */

constexpr int root = 0; /* root process rank */

/*
 * Timer variable(s).
 */
double elapsed, mintime, maxtime, avgtime;

int parse_cli(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    {
        /*
         * CombBLAS 2D processor grid uses MPI_COMM_WORLD
         * as its world communicator (comm).
         */
        Grid commgrid;

        comm = MPI_COMM_WORLD;
        commgrid.reset(new CommGrid(comm, 0, 0));
        myrank = commgrid->GetRank();
        nprocs = commgrid->GetSize();

        /*
         * Parse command-line from root process and broadcast.
         */
        if (parse_cli(argc, argv) != 0)
            goto err;

        MPI_Barrier(comm);
        elapsed = MPI_Wtime();

        FastaIndex index(fasta_fname, commgrid);
        Vector<String> myreads = index.GetMyReads();

        elapsed = MPI_Wtime() - elapsed;
        MPI_REDUCE(&elapsed, &mintime, 1, MPI_DOUBLE, MPI_MIN, root, comm);
        MPI_REDUCE(&elapsed, &maxtime, 1, MPI_DOUBLE, MPI_MAX, root, comm);
        MPI_REDUCE(&elapsed, &avgtime, 1, MPI_DOUBLE, MPI_SUM, root, comm);

        if (myrank == root)
        {
            avgtime /= nprocs;
            std::cout << "elapsed time: min=" << mintime << ", max=" << maxtime << ", avg=" << avgtime << " (seconds)" << std::endl;
        }

        MPI_Barrier(comm);

        goto done;
    }

err: returncode = -1;
done:
    MPI_Finalize();
    return returncode;
}

void usage(char const *prg)
{
    std::cerr << "Usage: " << prg << " [options] <reads.fa>\n"
              << "Options: -x INT   x-drop alignment threshold [" <<  xdrop_cutoff << "]\n"
              << "         -A INT   matching score ["             <<  mat          << "]\n"
              << "         -B INT   mismatch penalty ["           << -mis          << "]\n"
              << "         -G INT   gap penalty ["                << -gap          << "]\n"
              << "         -h       help message"
              << std::endl;
}

int parse_cli(int argc, char *argv[])
{
    int params[4] = {mat, mis, gap, xdrop_cutoff};
    int show_help = 0, fasta_provided = 1;

    if (myrank == root)
    {
        int c;

        while ((c = getopt(argc, argv, "x:A:B:G:h")) >= 0)
        {
            if      (c == 'A') params[0] =  atoi(optarg);
            else if (c == 'B') params[1] = -atoi(optarg);
            else if (c == 'G') params[2] = -atoi(optarg);
            else if (c == 'x') params[3] =  atoi(optarg);
            else if (c == 'h') show_help = 1;
        }
    }

    MPI_BCAST(params, 4, MPI_INT, root, comm);

    mat          = params[0];
    mis          = params[1];
    gap          = params[2];
    xdrop_cutoff = params[3];

    if (myrank == root && show_help)
        usage(argv[0]);

    MPI_BCAST(&show_help, 1, MPI_INT, root, comm);
    if (show_help) return -1;

    if (myrank == root && optind >= argc)
    {
        std::cerr << "error: missing FASTA file\n";
        usage(argv[0]);
        fasta_provided = 0;
    }

    MPI_BCAST(&fasta_provided, 1, MPI_INT, root, comm);
    if (!fasta_provided) return -1;

    int fnamelen;

    if (myrank == root)
    {
        fasta_fname = argv[optind];
        fnamelen = fasta_fname.size();

        std::cout << "-DKMER_SIZE="       << KMER_SIZE       << "\n"
                  << "-DLOWER_KMER_FREQ=" << LOWER_KMER_FREQ << "\n"
                  << "-DUPPER_KMER_FREQ=" << UPPER_KMER_FREQ << "\n"
        #ifdef USE_BLOOM
                  << "-DUSE_BLOOM\n"
        #endif
                  << "\n"
                  << "int mat = "          << mat          <<   ";\n"
                  << "int mis = "          << mis          <<   ";\n"
                  << "int gap = "          << gap          <<   ";\n"
                  << "int xdrop_cutoff = " << xdrop_cutoff <<   ";\n"
                  << "String fname = \""   << fasta_fname  << "\";\n"
                  << std::endl;
    }

    MPI_BCAST(&fnamelen, 1, MPI_INT, root, comm);

    if (myrank != root) fasta_fname.assign(fnamelen, '\0');

    MPI_BCAST(fasta_fname.data(), fnamelen, MPI_CHAR, root, comm);

    return 0;
}
