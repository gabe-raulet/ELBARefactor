#ifndef HYPER_LOG_LOG_H_
#define HYPER_LOG_LOG_H_

#include "common.h"
#include <cstdint>
#include <mpi.h>

class HyperLogLog
{
public:
    HyperLogLog(uint8_t bits);

    void Add(char const *s, size_t len);
    void Add(const String& s) { Add(s.c_str(), s.size()); }
    double Estimate() const;
    void Merge(const HyperLogLog& rhs);
    void ParallelMerge(MPI_Comm comm);

private:
    uint8_t bits; /* register bit width */
    uint32_t size; /* register size */
    double alpha_mm; /* alpha * m^2 */
    Vector<uint8_t> registers; /* registers */
};

#endif
