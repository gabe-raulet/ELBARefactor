#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "common.h"
#include "CommGrid.h"

typedef struct faidx_record { size_t len, pos, bases; } faidx_record_t;

class FastaIndex
{
public:
    FastaIndex(const String& fasta_fname, SharedPtr<CommGrid> commgrid);

    SharedPtr<CommGrid> getcommgrid() const { return commgrid; }
    const Vector<faidx_record_t>& getrecords() const { return records; }
    static Vector<String> GetMyReads(const FastaIndex& index);

    void PrintInfo() const;

    String GetFastaFilename() const { return fasta_fname; }
    String GetFaidxFilename() const { return fasta_fname + ".fai"; }

private:
    SharedPtr<CommGrid> commgrid;
    Vector<faidx_record_t> records;
    Vector<String> names;
    String fasta_fname;
};

#endif
