#ifndef LOGGER_H_
#define LOGGER_H_

#include "common.h"

void LogAll(const String mylog, SharedPtr<CommGrid> commgrid);
String ProcessorName(SharedPtr<CommGrid> commgrid);

class Logger
{
    std::unique_ptr<std::ostringstream> logstream, rootstream;
    SharedPtr<CommGrid> commgrid;
    int myrank, nprocs;
    MPI_Comm comm;

    String prefix();

public:
    Logger(SharedPtr<CommGrid> commgrid);
    void Flush(char const *label);
    void Flush(std::ostringstream& ss, int rank);
    std::ostringstream& operator()() { return *logstream; }
};


#endif
