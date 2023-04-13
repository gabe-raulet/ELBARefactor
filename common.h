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

template <class T>
using Vector = std::vector<T>;

template <class K, class V>
using Map = std::unordered_map<K, V>;

template <class T>
using Set = std::set<T>;

using String = std::string;

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

#endif
