#include "HashFuncs.h"

/* ----------- MURMURHASH ------------- */

std::uint32_t rotl32(std::uint32_t x, std::int8_t r)
{
    return (x << r) | (x >> (32 - r));
}

std::uint64_t rotl64(std::uint64_t x, std::int8_t r)
{
    return (x << r) | (x >> (64 - r));
}

#define ROTL32(x,y) rotl32(x,y) /* GRGR: why is it necessary to use these macros instead of the original functions? */
#define ROTL64(x,y) rotl64(x,y)

#define BIG_CONSTANT(x) (x##ULL)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

#define getblock(p, i) (p[i])

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

std::uint64_t fmix64(std::uint64_t k)
{
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;

    return k;
}

void murmurhash3_x64_128(const void *key, const std::uint32_t len, const std::uint32_t seed, void *hashval)
{
    const std::uint8_t *data = (const std::uint8_t*)key;
    const std::uint32_t nblocks = len / 16;
    std::int32_t i;

    std::uint64_t h1 = seed;
    std::uint64_t h2 = seed;

    std::uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
    std::uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

    const std::uint64_t *blocks = (const std::uint64_t *)(data);

    for(i = 0; i < nblocks; i++)
    {
        std::uint64_t k1 = getblock(blocks,i*2+0);
        std::uint64_t k2 = getblock(blocks,i*2+1);

        k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

        h1 = ROTL64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

        k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

        h2 = ROTL64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
    }

    const std::uint8_t *tail = (const std::uint8_t*)(data + nblocks*16);

    std::uint64_t k1 = 0;
    std::uint64_t k2 = 0;

    switch(len & 15)
    {
    case 15: k2 ^= (std::uint64_t)(tail[14]) << 48;
    case 14: k2 ^= (std::uint64_t)(tail[13]) << 40;
    case 13: k2 ^= (std::uint64_t)(tail[12]) << 32;
    case 12: k2 ^= (std::uint64_t)(tail[11]) << 24;
    case 11: k2 ^= (std::uint64_t)(tail[10]) << 16;
    case 10: k2 ^= (std::uint64_t)(tail[ 9]) << 8;
    case  9: k2 ^= (std::uint64_t)(tail[ 8]) << 0;
        k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= (std::uint64_t)(tail[ 7]) << 56;
    case  7: k1 ^= (std::uint64_t)(tail[ 6]) << 48;
    case  6: k1 ^= (std::uint64_t)(tail[ 5]) << 40;
    case  5: k1 ^= (std::uint64_t)(tail[ 4]) << 32;
    case  4: k1 ^= (std::uint64_t)(tail[ 3]) << 24;
    case  3: k1 ^= (std::uint64_t)(tail[ 2]) << 16;
    case  2: k1 ^= (std::uint64_t)(tail[ 1]) << 8;
    case  1: k1 ^= (std::uint64_t)(tail[ 0]) << 0;
        k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
    };

    h1 ^= len; h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    ((std::uint64_t*)hashval)[0] = h1;
    ((std::uint64_t*)hashval)[1] = h2;
}

/* ------------ THOMAS WANG INVERTIBLE INTEGER HASH FUNCTION ------------ */

void wang_hash_64bits(const void *key, void *hashval)
{
    std::uint64_t _key = *((std::uint64_t const *)key);

    _key = (~_key) + (_key << 21);
    _key = _key ^ (_key >> 24);
    _key = (_key + (_key << 3)) + (_key << 8);
    _key = _key ^ (_key >> 14);
    _key = (_key + (_key << 2)) + (_key << 4);
    _key = _key ^ (_key >> 28);
    _key = _key + (_key << 31);

    *((std::uint64_t*)hashval) = _key;
}

/* reference: https://naml.us/post/inverse-of-a-hash-function/ */

void wang_inverse_hash_64bits(const void *hashval, void *key)
{
    std::uint64_t tmp;
    std::uint64_t _hashval = *((std::uint64_t const *)hashval);

    tmp = _hashval - (_hashval << 31);
    _hashval = _hashval - (tmp << 31);

    tmp = _hashval ^ _hashval >> 28;
    _hashval = _hashval ^ tmp >> 28;

    _hashval *= BIG_CONSTANT(14933078535860113213);

    tmp = _hashval ^ _hashval >> 14;
    tmp = _hashval ^ tmp >> 14;
    tmp = _hashval ^ tmp >> 14;
    _hashval = _hashval ^ tmp >> 14;

    _hashval *= BIG_CONSTANT(15244667743933553977);

    tmp = _hashval ^ _hashval >> 24;
    _hashval = _hashval ^ tmp >> 24;

    tmp = ~_hashval;
    tmp = ~(_hashval - (tmp << 21));
    tmp = ~(_hashval - (tmp << 21));
    _hashval = ~(_hashval - (tmp << 21));

    *((std::uint64_t*)key) = _hashval;
}

/* ------------ MURMURHASH WRAPPERS -------------- */

void murmurhash3_128bits(const void *key, std::uint32_t numbytes, void *hashval)
{
    murmurhash3_x64_128(key, numbytes, 313, hashval);
}

void murmurhash3_64bits(const void *key, std::uint32_t numbytes, void *hashval)
{
    std::uint64_t vals[2];
    murmurhash3_x64_128(key, numbytes, 313, (void*)vals);
    *((std::uint64_t*)hashval) = vals[0];
}

void murmurhash3_32bits(const void *key, std::uint32_t numbytes, void *hashval)
{
    std::uint64_t vals[2];
    murmurhash3_x64_128(key, numbytes, 313, (void*)vals);
    *((std::uint32_t*)hashval) = vals[0] & ((1ULL<<32)-1);
}

std::uint32_t murmurhash3(const void *key, std::size_t len, std::uint32_t seed)
{
    std::uint64_t vals[2];
    murmurhash3_x64_128(key, len, seed, (void*)vals);
    return (std::uint32_t)(vals[0] & ((1ULL<<32)-1));
}

