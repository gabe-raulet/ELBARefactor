#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "common.h"

typedef struct faidx_record { size_t len, pos, bases; } faidx_record_t;

class FastaIndex
{
public:
    FastaIndex(const String& fasta_fname, SharedPtr<CommGrid> commgrid);

    SharedPtr<CommGrid> getcommgrid() const { return commgrid; }
    Vector<String> GetAllReads() { return GetReadsFromRecords(allrecords); }
    Vector<String> GetMyReads();

    String GetFastaFilename() const { return fasta_fname; }
    String GetFaidxFilename() const { return fasta_fname + ".fai"; }

private:
    SharedPtr<CommGrid> commgrid;
    Vector<faidx_record_t> myrecords, allrecords;
    Vector<String> names;
    String fasta_fname;

    Vector<String> GetReadsFromRecords(const Vector<faidx_record_t>& records);
};

#endif
