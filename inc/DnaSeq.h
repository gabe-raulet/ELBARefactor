#ifndef DnaSeq_H_
#define DnaSeq_H_

#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

class DnaSeq
{
public:
    DnaSeq() = default;
    DnaSeq(char const *sequence, size_t len);
    DnaSeq(char const *sequence) : DnaSeq(sequence, strlen(sequence)) {}
    DnaSeq(std::string const& sequence) : DnaSeq(sequence.c_str(), sequence.size()) {}
    DnaSeq(const DnaSeq& rhs) : words(rhs.words), remain(rhs.remain) {}

    DnaSeq& operator=(DnaSeq rhs)
    {
        std::swap(words, rhs.words);
        std::swap(remain, rhs.remain);
        return *this;
    }

    std::string ascii() const;
    DnaSeq substr(size_t pos, size_t count) const;
    size_t size() const { return words.size() * 32 - remain; }
    size_t numbytes() const { return words.size() * 8; }

    int operator[](size_t i) const;
    bool operator==(const DnaSeq& rhs);
    bool operator<(const DnaSeq& rhs);
    bool operator!=(const DnaSeq& rhs) { return !(*this == rhs); }

    static char    getcodechar(int c)  { return chartab[c]; }
    static uint8_t getcharcode(char c) { return codetab[(int)c]; }

    friend std::ostream& operator<<(std::ostream& stream, const DnaSeq& s)
    {
        stream << s.ascii();
        return stream;
    }

private:
    std::vector<uint64_t> words; /* each word encodes up to 32 nucleotides */
    int remain;                  /* number of nucleotides in last word */

    static constexpr char chartab[4] = {'A', 'C', 'G', 'T'};
    static constexpr uint8_t codetab[256] =
    {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
};

#endif
