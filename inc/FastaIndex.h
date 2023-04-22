#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "common.h"

class FastaIndex
{
public:
    typedef struct faidx_record { size_t len, pos, bases; } faidx_record_t;

    FastaIndex(const String& fasta_fname, Grid commgrid);

    Vector<String> GetAllReads();
    Vector<String> GetMyReads();

    String GetFastaFilename() const { return fasta_fname; }
    String GetFaidxFilename() const { return fasta_fname + ".fai"; }

private:
    Grid commgrid;
    Vector<faidx_record_t> myrecords, allrecords;
    Vector<String> names;
    String fasta_fname;

    Vector<String> GetReadsFromRecords(const Vector<faidx_record_t>& records);
};

#endif
