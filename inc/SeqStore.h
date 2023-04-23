#ifndef SEQ_STORE_H_
#define SEQ_STORE_H_

#include "common.h"
#include "FastaIndex.h"

class SeqStore
{
public:
    SeqStore(const FastaIndex& index);

private:
    Grid commgrid;
    Vector<char> buf;
    Vector<size_t> offsets;
    Vector<size_t> readlens;
};

#endif
