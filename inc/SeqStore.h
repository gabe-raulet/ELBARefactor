#ifndef SEQ_STORE_H_
#define SEQ_STORE_H_

#include "common.h"
#include "FastaIndex.h"

class SeqStore
{
public:
    SeqStore(const FastaIndex& index);
    ~SeqStore();

    char const* GetSequence(size_t id) const { return buf + offsets[id]; }
    size_t GetMyNumReads() const { return offsets.size(); }
    size_t GetNumReads() const;
    Vector<size_t> GetGlobalReadCounts() const;
    Vector<size_t> GetGlobalReadOffsets() const;

private:
    Grid commgrid;
    char *buf;
    Vector<size_t> offsets;
    Vector<size_t> readlens;
};

#endif
