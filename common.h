#ifndef COMMON_H_
#define COMMON_H_

#include <vector>
#include <array>
#include <string>
#include <tuple>
#include <unordered_map>
#include <set>
#include <memory>
#include <cstdint>
#include <mpi.h>

using String = std::string;

template <class T>
using Vector = std::vector<T>;

template <class K, class V>
using Map = std::unordered_map<K, V>;

template <class T>
using Set = std::set<T>;

template <class T>
using SharedPtr = std::shared_ptr<T>;

template <class T, std::size_t N>
using Array = std::array<T, N>;

template <class T1, class T2>
using Pair = std::pair<T1, T2>;

template <class... Types>
using Tuple = std::tuple<Types...>;

template <class T>
using Hash = std::hash<T>;

// #if MPI_VERSION == 3
// #define MPI_HAS_LARGE_COUNTS 0

// typedef int MPI_Count_type;
// typedef int MPI_Displ_type;

// #define MPI_COUNT_TYPE MPI_INT

// #define MPI_ALLTOALL MPI_Alltoall
// #define MPI_ALLTOALLV MPI_Alltoallv

// #elif MPI_VERSION == 4
// #define MPI_HAS_LARGE_COUNTS 1

// typedef MPI_Count MPI_Count_type;
// typedef MPI_Aint MPI_Displ_type;

// #define MPI_COUNT_TYPE MPI_COUNT

// #define MPI_ALLTOALL MPI_Alltoall_c
// #define MPI_ALLTOALLV MPI_Alltoallv_c

// #else
// #error "MPI version should either be 3 or 4."
// #endif

#endif
