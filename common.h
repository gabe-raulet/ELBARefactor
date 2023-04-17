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
#include <limits>
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

#if MPI_VERSION == 3
#define MPI_HAS_LARGE_COUNTS 0

typedef int MPI_Count_type;
typedef int MPI_Displ_type;

#define MPI_COUNT_TYPE MPI_INT

#define MPI_ALLTOALL MPI_Alltoall
#define MPI_ALLTOALLV MPI_Alltoallv
#define MPI_SCATTER MPI_Scatter
#define MPI_SCATTERV MPI_Scatterv
#define MPI_GATHER MPI_Gather
#define MPI_GATHERV MPI_Gatherv
#define MPI_FILE_READ_AT_ALL MPI_File_read_at_all

#elif MPI_VERSION == 4
#define MPI_HAS_LARGE_COUNTS 1

typedef MPI_Count MPI_Count_type;
typedef MPI_Aint MPI_Displ_type;

#define MPI_COUNT_TYPE MPI_COUNT

#define MPI_ALLTOALL MPI_Alltoall_c
#define MPI_ALLTOALLV MPI_Alltoallv_c
#define MPI_SCATTER MPI_Scatter_c
#define MPI_SCATTERV MPI_Scatterv_c
#define MPI_GATHER MPI_Gather_c
#define MPI_GATHERV MPI_Gatherv_c
#define MPI_FILE_READ_AT_ALL MPI_File_read_at_all_c

#else
#error "MPI version should either be 3 or 4."
#endif

#ifndef MPI_SIZE_T
static_assert(std::numeric_limits<size_t>::max() == std::numeric_limits<unsigned long>::max());
#define MPI_SIZE_T MPI_UNSIGNED_LONG
#endif

#endif
