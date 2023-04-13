#ifndef HASHFUNCS_H_
#define HASHFUNCS_H_

#include <cstdint>

void wang_hash_64bits(const void *key, void *hashval);
void wang_inverse_hash_64bits(const void *hashval, void *key);

void murmurhash3_128bits(const void *key, std::uint32_t numbytes, void *hashval);
void murmurhash3_64bits(const void *key, std::uint32_t numbytes, void *hashval);
void murmurhash3_32bits(const void *key, std::uint32_t numbytes, void *hashval);

std::uint32_t murmurhash3(const void *key, std::size_t len, std::uint32_t seed);

#endif
