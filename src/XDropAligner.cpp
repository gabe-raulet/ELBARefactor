#include "XDropAligner.h"
#include <algorithm>

XDropAligner::XDropAligner(const Vector<String>& allreads, int mat, int mis, int gap, int dropoff, SharedPtr<CommGrid> commgrid) : allreads(allreads), mat(mat), mis(mis), gap(gap), dropoff(dropoff), commgrid(commgrid) {}

void XDropAligner::Apply(CT<ReadOverlap>::PSpParMat& B)
{
    CT<ReadOverlap>::PSpDCCols *spSeq = B.seqptr();

    uint64_t total_rows = B.getnrow();
    uint64_t total_cols = B.getncol();
    int procrows = commgrid->GetGridRows();
    int proccols = commgrid->GetGridCols();
    uint64_t myrowoffset = commgrid->GetRankInProcCol() * (total_rows / procrows);
    uint64_t mycoloffset = commgrid->GetRankInProcRow() * (total_cols / proccols);

    std::ostringstream ss;

    for (auto colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
    {
        for (auto nzit = spSeq->begnz(colit); nzit != spSeq->endnz(colit); ++nzit)
        {
            uint64_t row = nzit.rowid() + myrowoffset;
            uint64_t col = colit.colid() + mycoloffset;
            ReadOverlap& overlap = nzit.value();
        }
    }

    for (int i = 0; i < commgrid->GetSize(); ++i)
    {
        if (commgrid->GetRank() == i)
        {
            std::cout << ss.str() << std::flush;
        }

        MPI_Barrier(commgrid->GetWorld());
    }
}

char revcomp(char c)
{
    switch (c)
    {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'T';
        default: break;
    }
    return 'N';
}



String GetReverseComplement(const String& s)
{
    String sequence(s);
    std::reverse(sequence.begin(), sequence.end());
    std::transform(sequence.cbegin(), sequence.cend(), sequence.begin(), revcomp);
    return sequence;
}

int XDropAligner::ExtendLeft(const String& seqQ, const String& seqT, const XSeed& xseed, int& begQ_ext, int& begT_ext)
{
    XSeed result(xseed);

    int lscore = ExtendSeedOneDirection(seqQ, seqT, true, result);

    begQ_ext = result.begQ;
    begT_ext = result.begT;

    return lscore;
}

int XDropAligner::ExtendRight(const String& seqQ, const String& seqT, const XSeed& xseed, int& endQ_ext, int& endT_ext)
{
    XSeed result(xseed);

    int rscore = ExtendSeedOneDirection(seqQ, seqT, false, result);

    endQ_ext = result.endQ;
    endT_ext = result.endT;

    return rscore;
}

XDropAligner::XSeed XDropAligner::SeedAndExtend(const String& seqQ, const String& seqT, int begQ, int begT)
{
    const String seqTr = GetReverseComplement(seqT);
    XSeed xseed, result;

    int lenQ = seqQ.size();
    int lenT = seqT.size();

    assert(begQ >= 0 && begQ + KMER_SIZE <= lenQ);
    assert(begT >= 0 && begT + KMER_SIZE <= lenT);

    bool rc = (seqQ[begQ + (KMER_SIZE>>1)] != seqT[begT + (KMER_SIZE>>1)]);

    xseed.begQ = begQ;
    xseed.endQ = xseed.begQ + KMER_SIZE;

    xseed.begT = rc? lenT - begT - KMER_SIZE : begT;
    xseed.endT = xseed.begT + KMER_SIZE;

    xseed.rc = rc;

    int begQ_ext, begT_ext, lscore;
    int endQ_ext, endT_ext, rscore;

    lscore = ExtendLeft(seqQ, seqT, xseed, begQ_ext, begT_ext);
    rscore = ExtendRight(seqQ, seqT, xseed, endQ_ext, endT_ext);

    int score = lscore + rscore + mat * KMER_SIZE;

    result.begQ = begQ_ext;
    result.endQ = endQ_ext;

    result.begT = rc? lenT - endT_ext : begT_ext;
    result.endT = rc? lenT - begT_ext : endT_ext;

    result.rc = rc;
    result.score = score;

    return result;
}

int XDropAligner::ExtendSeedOneDirection(const String& seqQ, const String& seqT, bool extleft, XSeed& xseed)
{
    int lenQ = seqQ.size();
    int lenT = seqT.size();

    const String seqTr = GetReverseComplement(seqT);

    int lenQ_ext = extleft? xseed.begQ : lenQ - xseed.endQ;
    int lenT_ext = extleft? xseed.begT : lenT - xseed.endT;

    int cols = lenQ_ext + 1;
    int rows = lenT_ext + 1;

    if (rows == 1 || cols == 1) return 0;

    constexpr int int_min = std::numeric_limits<int>::min();

    int len = 2 * std::max(cols, rows);
    int min_err_score = int_min / len;
    gap = std::max(gap, min_err_score);
    mis = std::max(mis, min_err_score);
    int undef = int_min - gap;

    std::vector<int> ad1, ad2, ad3, tmp;

    int min_col = 1, max_col = 2;
    int offset1 = 0, offset2 = 0, offset3 = 0;

    ad2.push_back(0);
    int best_ext_col = 0, best_ext_row = 0, best_ext_score = 0;

    ad3.resize(2);
    ad3[0] = ad3[1] = (-gap > dropoff)? undef : gap;

    int ad_no = 1, best = 0;
    int offsetQ = xseed.endQ;
    int offsetT = xseed.endT;

    while (min_col < max_col)
    {
        ++ad_no;
        tmp = std::move(ad1);
        ad1 = std::move(ad2);
        ad2 = std::move(ad3);
        ad3 = std::move(tmp);

        offset1 = offset2;
        offset2 = offset3;
        offset3 = min_col - 1;

        ad3.resize(max_col + 1 - offset3);
        ad3[0] = ad3[max_col - offset3] = undef;

        if (ad_no * gap > best - dropoff)
        {
            if (offset3 == 0) ad3[0] = ad_no * gap;
            if (ad_no - max_col == 0) ad3[max_col - offset3] = ad_no * gap;
        }

        int ad_best = ad_no * gap;

        for (int col = min_col; col < max_col; ++col)
        {
            int i3 = col - offset3;
            int i2 = col - offset2;
            int i1 = col - offset1;

            int posQ, posT;

            posQ = extleft? cols - 1 - col : col - 1 + offsetQ;
            posT = extleft? rows - 1 + col - ad_no : ad_no - col - 1 + offsetT;

            int temp = std::max(ad2[i2-1], ad2[i2]) + gap;
            int temp2 = ad1[i1-1] + ((seqQ[posQ] == (xseed.rc? seqTr[posT] : seqT[posT]))? mat : mis);
            temp = std::max(temp, temp2);

            if (temp < best - dropoff)
            {
                ad3[i3] = undef;
            }
            else
            {
                ad3[i3] = temp;
                ad_best = std::max(ad_best, temp);
            }

            if (temp > best)
            {
                best_ext_col = col;
                best_ext_row = ad_no - best_ext_col;
                best_ext_score = ad3[best_ext_col - offset3];
                assert(best_ext_score == temp);
            }
        }

        best = std::max(best, ad_best);

        while (min_col - offset3 < ad3.size() && ad3[min_col - offset3] == undef &&
               min_col - offset2 - 1 < ad2.size() && ad2[min_col - offset2 - 1] == undef)
        {
            ++min_col;
        }

        while (max_col - offset3 > 0 && ad3[max_col - offset3 - 1] == undef && ad2[max_col - offset2 - 1] == undef)
            --max_col;

        ++max_col;

        min_col = std::max(min_col, ad_no + 2 - rows);
        max_col = std::min(max_col, cols);
    }

    int ext_col = ad3.size() + offset3 - 2;
    int ext_row = ad_no - ext_col;
    int ext_score = ad3[ext_col - offset3];

    if (ext_score == undef)
    {
        if (ad2[ad2.size()-2] != undef)
        {
            ext_col = ad2.size() + offset2 - 2;
            ext_row = ad_no - 1 - ext_col;
            ext_score = ad2[ext_col - offset2];
        }
        else if (ad2.size() > 2 && ad2[ad2.size()-3] != undef)
        {
            ext_col = ad2.size() + offset2 - 3;
            ext_row = ad_no - 1 - ext_col;
            ext_score = ad2[ext_col - offset3];
        }
    }

    if (ext_score == undef)
    {
        for (int i = 0; i < ad1.size(); ++i)
        {
            if (ad1[i] > ext_score)
            {
                ext_score = ad1[i];
                ext_col = i + offset1;
                ext_row = ad_no - 2 - ext_col;
            }
        }
    }

    if (best_ext_score != undef)
    {
        if (extleft)
        {
            xseed.begT -= best_ext_row;
            xseed.begQ -= best_ext_col;
        }
        else
        {
            xseed.endT += best_ext_row;
            xseed.endQ += best_ext_col;
        }
    }

    return best_ext_score;

}

void XDropAligner::ClassifyAlignment(const XSeed& ai, int lenQ, int lenT, OverlapClass& kind)
{
    if (ai.score <= 0)
    {
        kind = BAD_ALIGNMENT;
        return;
    }

    int begTr = ai.rc? lenT - ai.endT : ai.begT;
    int endTr = ai.rc? lenT - ai.begT : ai.endT;

    int maplen = ((ai.endT - ai.begT) + (ai.endQ - ai.begQ)) / 2;
    int overhang = std::min(ai.begQ, begTr) + std::min(lenQ - ai.endQ, lenT - endTr);
    int overlap = maplen + overhang;

    float my_thr = 0.9 * (0.99 * (overlap + 0.0));

    if (ai.begQ <= begTr && lenQ - ai.endQ <= lenT - endTr)
    {
        kind = FIRST_CONTAINED;
    }
    else if (ai.begQ >= begTr && lenQ - ai.endQ >= lenT - endTr)
    {
        kind = SECOND_CONTAINED;
    }
    else if (ai.score < my_thr || overlap < 500)
    {
        kind = BAD_ALIGNMENT;
    }
    else if (ai.begQ > begTr)
    {
        kind = FIRST_TO_SECOND_OVERLAP;
    }
    else
    {
        kind = SECOND_TO_FIRST_OVERLAP;
    }
}
