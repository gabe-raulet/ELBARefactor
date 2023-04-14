#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cassert>
#include "common.h"
#include "Kmer.h"
#include "KmerComm.h"
#include "CommGrid.h"
#include "FastaIndex.h"

int main(int argc, char *argv[])
{
    TKmer::SetKmerSize(51);

    String myread = "ACTACTGGCAAAAAAAAAATGCATCGATCGATCGCGCGTATGGTTTTTTTTTCATCGCGGCCCCCCCGGGCTAGCATCGACTAGCTGGCGGGATTCGGGCGCGCA";

    Vector<TKmer> realkmers = TKmer::GetKmers(myread);
    Vector<TKmer> repkmers = TKmer::GetRepKmers(myread);

    assert(realkmers.size() == repkmers.size());

    for (size_t i = 0; i < realkmers.size(); ++i)
    {
        std::cout << realkmers[i] << "\t" << repkmers[i] << std::endl;
        assert(repkmers[i] < realkmers[i] || repkmers[i] == realkmers[i]);
        assert(repkmers[i] == realkmers[i].GetRep());
    }

    return 0;
}
