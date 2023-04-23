#include "DnaSeq.h"
#include <cassert>
#include <cstring>

DnaSeq::DnaSeq(char const *sequence, size_t len)
{
    size_t numwords = (len + 31) / 32;
    remain = (numwords * 32) - len;
    assert(remain < 32);

    words.resize(numwords);

    uint64_t word;
    size_t w = 0;
    char const *sw = sequence;

    while (w < numwords)
    {
        word = 0;
        int left = (w != numwords-1? 32 : 32-remain);

        for (int i = 0; i < left; ++i)
        {
            uint8_t code = DnaSeq::getcharcode(sw[i]);
            uint64_t shift = static_cast<uint64_t>(code) << (62 - (2*i));
            word |= shift;
        }

        words[w++] = word;
        sw += 32;
    }
}

std::string DnaSeq::ascii() const
{
    size_t len = size();
    uint64_t const *wb = words.data();
    std::vector<char> s(len);

    for (size_t i = 0; i < len; ++i)
    {
        int code = (*wb >> (62 - (2*(i % 32)))) & 3;
        s[i] = DnaSeq::getcodechar(code);

        if ((i+1) % 32 == 0)
            wb++;
    }

    return std::string(s.begin(), s.end());
}

DnaSeq DnaSeq::substr(size_t pos, size_t count) const
{
    return DnaSeq(ascii().substr(pos, count));
}

int DnaSeq::operator[](size_t i) const
{
    uint64_t word = words[i/32];
    int shift = 62 - (2 * (i%32));
    int code = (word >> shift)&3;
    return code;
}

bool DnaSeq::operator==(const DnaSeq& rhs)
{
    if (size() != rhs.size())
        return false;

    for (size_t i = 0; i < rhs.size(); ++i)
        if ((*this)[i] != rhs[i])
            return false;

    return true;
}

bool DnaSeq::operator<(const DnaSeq& rhs)
{
    size_t len = std::min(size(), rhs.size());
    size_t i;

    for (i = 0; (*this)[i] == rhs[i] && i < len; ++i);

    if (i == len) return false;

    return ((*this)[i] < rhs[i]);
}
