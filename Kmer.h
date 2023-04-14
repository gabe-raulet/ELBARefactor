#ifndef KMER_H_
#define KMER_H_

#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdint>
#include <cstring>
#include "common.h"
#include "HyperLogLog.h"

template <int N_LONGS>
class Kmer
{
public:

    static constexpr int N_BYTES = 8 * N_LONGS;
    static inline int k = 32*N_LONGS - 1;

    typedef Array<uint64_t, N_LONGS> MERARR;
    typedef Array<uint8_t,  N_BYTES> BYTEARR;

    static void SetKmerSize(int kmer_size)
    {
        assert((kmer_size & 1) && ((32*(N_LONGS-1)) <= kmer_size) && ((kmer_size < 32*N_LONGS)));
        k = kmer_size;
    }

    Kmer();
    Kmer(char const *s);
    Kmer(const Kmer& o);

    Kmer& operator=(const Kmer& o);

    String GetString() const;

    bool operator<(const Kmer& o) const;
    bool operator>(const Kmer& o) const;
    bool operator==(const Kmer& o) const;
    bool operator!=(const Kmer& o) const;

    Kmer GetExtension(char const nt) const;
    Kmer GetTwin() const;
    Kmer GetRep() const;

    uint64_t GetHash() const;
    const void* GetBytes() const { return reinterpret_cast<const void*>(longs.data()); }

    int GetInferredOwner(int nprocs) const;

    void CopyDataInto(void *mem) const { std::memcpy(mem, longs.data(), N_BYTES); }
    void CopyDataFrom(const void *mem) { std::memcpy(longs.data(), mem, N_BYTES); }

    static Vector<Kmer> GetKmers(const String& s);
    static Vector<Kmer> GetRepKmers(const String& s);

    static void InsertIntoHLL(const String& s, HyperLogLog& hll);

    template <int N>
    friend std::ostream& operator<<(std::ostream& os, const Kmer<N>& kmer);

private:

    union { MERARR  longs;
            BYTEARR bytes; };

    void set_kmer(char const *s, bool const revcomp = false);
};

template <int N_LONGS>
std::ostream& operator<<(std::ostream& os, const Kmer<N_LONGS>& kmer)
{
    os << kmer.k << "-mer(" << kmer.GetString() << ")";
    return os;
}

namespace std
{
    template <int N_LONGS> struct hash<Kmer<N_LONGS>>
    {
        size_t operator()(const Kmer<N_LONGS>& kmer) const
        {
            auto myhash = kmer.GetHash();
            return myhash;
        }
    };

    template <int N_LONGS> struct less<Kmer<N_LONGS>>
    {
        bool operator()(const Kmer<N_LONGS>& k1, const Kmer<N_LONGS>& k2) const
        {
            return k1 < k2;
        }
    };
}

#include "Kmer.cpp"

typedef Kmer<1> TKmer;

#endif
