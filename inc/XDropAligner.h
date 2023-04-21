#ifndef XDROP_ALIGNMENT_H_
#define XDROP_ALIGNMENT_H_

#include "common.h"
#include "ReadOverlap.h"

class XDropAligner
{
public:
    XDropAligner(const Vector<String>& allreads, int mat, int mis, int gap, int dropoff, SharedPtr<CommGrid> commgrid);
    void Apply(CT<ReadOverlap>::PSpParMat& B);

    enum OverlapClass
    {
        BAD_ALIGNMENT, /* bad alignment */
        FIRST_CONTAINED,
        SECOND_CONTAINED,
        FIRST_TO_SECOND_OVERLAP,
        SECOND_TO_FIRST_OVERLAP
    };

    struct XSeed
    {
        int begQ, endQ, begT, endT, score;
        bool rc;

        XSeed() : begQ(0), endQ(0), begT(0), endT(0), score(-1), rc(false) {}

        XSeed(const XSeed& rhs) = default;
    };

    static void ClassifyAlignment(const XSeed& ai, int lenQ, int lenT, OverlapClass& kind);
    XSeed SeedAndExtend(const String& seqQ, const String& seqT, int begQ, int begT);

private:
    const Vector<String>& allreads;
    int mat, mis, gap, dropoff;
    SharedPtr<CommGrid> commgrid;

    int ExtendSeedOneDirection(const String& seqQ, const String& seqT, bool extleft, XSeed& xseed);
    int ExtendLeft(const String& seqQ, const String& seqT, const XSeed& xseed, int& begQ_ext, int& begT_ext);
    int ExtendRight(const String& seqQ, const String& seqT, const XSeed& xseed, int& endQ_ext, int& endT_ext);
};


#endif
