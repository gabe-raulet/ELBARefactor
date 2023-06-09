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

#ifndef KMER_SIZE
#error "KMER_SIZE must be defined"
#endif

#if (KMER_SIZE <= 0) || (KMER_SIZE >= 96) || (!(KMER_SIZE&1))
#error "KMER_SIZE must be in the range (0,96) and must be odd"
#endif

template <int N_LONGS>
class Kmer
{
public:

    static_assert(N_LONGS != 0);

    static constexpr int N_BYTES = 8 * N_LONGS;

    typedef Array<uint64_t, N_LONGS> MERARR;
    typedef Array<uint8_t,  N_BYTES> BYTEARR;

    Kmer();
    Kmer(char const *s);
    Kmer(const void *mem);
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

    void CopyDataInto(void *mem) const { std::memcpy(mem, longs.data(), N_BYTES); }
    void CopyDataFrom(const void *mem) { std::memcpy(longs.data(), mem, N_BYTES); }

    static Vector<Kmer> GetKmers(const String& s);
    static Vector<Kmer> GetRepKmers(const String& s);

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
    os << KMER_SIZE << "-mer(" << kmer.GetString() << ")";
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

using TKmer = typename std::conditional<(KMER_SIZE <= 32), Kmer<1>,
              typename std::conditional<(KMER_SIZE <= 64), Kmer<2>,
              typename std::conditional<(KMER_SIZE <= 96), Kmer<3>, Kmer<0>>::type>::type>::type;

#endif
