#ifndef FASTA_DATA_H_
#define FASTA_DATA_H_

#include "common.h"
#include "FastaIndex.h"
#include <cassert>

class FastaData
{
public:
    FastaData(Grid commgrid) : commgrid(commgrid) {}
    FastaData(const FastaIndex& index);
    FastaData(const FastaData& rhs);
    size_t numreads() const { assert(offsets.size() == readlens.size()); return offsets.size(); }
    size_t getfirstid() const { return firstid; }
    FastaData getrange(size_t pos, size_t count) const;
    FastaData& operator+=(const FastaData& rhs);
    String getsequence(size_t localid) const;
    Grid getcommgrid() const { return commgrid; }

    ~FastaData() { delete[] buf; }

private:
    Grid commgrid;
    uint8_t *buf; /* buffer storing local read sequences (compressed) */
    Vector<size_t> offsets; /* offsets within compressed buf of first byte in a read sequence */
    Vector<size_t> readlens; /* uncompressed read lengths */
    size_t firstid; /* these represent global read ids */
    size_t bufsize;
};

class DistributedFastaData
{
public:
    DistributedFastaData(const FastaData& localdata);

    FastaData CollectRowReads() const;
    FastaData CollectColReads() const;

private:
    Grid commgrid;
    FastaData localdata;
    Vector<size_t> idoffsets;
};


#endif
