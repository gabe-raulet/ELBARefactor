#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "common.h"
#include "DnaSeq.h"

class FastaIndex
{
public:
    typedef struct faidx_record { size_t len, pos, bases; } faidx_record_t;

    FastaIndex(const String& fasta_fname, Grid commgrid, bool membalanced = false);

    Grid getcommgrid() const { return commgrid; }

    size_t GetNumRecords() const { return records.size(); }

    const Vector<faidx_record_t>& getmyrecords() const { return myrecords; }

    Vector<DnaSeq> GetMyReads();

    String GetFastaFilename() const { return fasta_fname; }
    String GetFaidxFilename() const { return fasta_fname + ".fai"; }

private:
    Grid commgrid;
    Vector<faidx_record_t> myrecords, records;
    Vector<String> names;
    String fasta_fname;

    Vector<DnaSeq> GetReadsFromRecords(const Vector<faidx_record_t>& records);

    MPI_Count_type get_idbalanced_partition(Vector<MPI_Count_type>& sendcounts);
    MPI_Count_type get_membalanced_partition(Vector<MPI_Count_type>& sendcounts);
};

#endif
