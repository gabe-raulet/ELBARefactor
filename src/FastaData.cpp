#include "FastaData.h"
#include "Logger.h"

FastaData::FastaData(const FastaIndex& index) : commgrid(index.getcommgrid())
{
    /*
     * Will skip convention of adding 'my' to the front of variables referring
     * to local processor information since that is implied by the nature
     * of this class' purpose.
     */
    size_t numreads;
    char *readbuf;
    MPI_Offset startpos, endpos, filesize, readbufsize;
    MPI_File fh;

    const auto& records = index.getmyrecords();
    numreads = records.size();
    startpos = records.front().pos;
    endpos = records.back().pos + records.back().len + (records.back().len / records.back().bases);
    MPI_File_open(commgrid->GetWorld(), index.GetFastaFilename().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_get_size(fh, &filesize);
    if (endpos > filesize) endpos = filesize;
    readbufsize = endpos - startpos;
    readbuf = new char[readbufsize];
    MPI_FILE_READ_AT_ALL(fh, startpos, readbuf, readbufsize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    /*
     * For a sequence of length l, floor((l+3)/4) bytes are needed to encode it. Let
     * L = l[0] + l[1] + ... + l[N-1], where l[i] is the length of the ith sequence, N is
     * the total number of sequences, and hence L is the total sum of all the sequence
     * lengths. Then the total number of bytes needed for the write buffer is
     *
     * Sum{0 <= i <= N-1}[floor((l[i]+3)/4)] <= (1/4) * Sum{0 <= i <= N-1}[l[i] + 4]
     *                                       <= (1/4) * (L + 4N)
     *                                        = (L/4) + N
     */

    readlens.reserve(numreads);
    offsets.reserve(numreads);

    size_t totbases = 0, maxlen = 0;
    for (auto itr = records.cbegin(); itr != records.cend(); ++itr)
    {
        totbases += itr->len;
        maxlen = std::max(itr->len, maxlen);
        readlens.push_back(itr->len);
    }

    bufsize = static_cast<size_t>(std::ceil((totbases/4.0) + numreads));
    buf = new uint8_t[bufsize];

    char *tmpbuf = new char[maxlen];

    size_t byteoffset = 0;

    MPI_Barrier(commgrid->GetWorld());
    double elapsed = MPI_Wtime();

    for (auto itr = records.cbegin(); itr != records.cend(); ++itr)
    {
        size_t locpos = 0;
        ptrdiff_t chunkpos = itr->pos - startpos;
        ptrdiff_t remain = itr->len;
        char *writeptr = tmpbuf;

        while (remain > 0)
        {
            size_t cnt = std::min(itr->bases, static_cast<size_t>(remain));
            std::memcpy(writeptr, &readbuf[chunkpos + locpos], cnt);
            writeptr += cnt;
            remain -= cnt;
            locpos += (cnt+1);
        }

        size_t numbytes = (itr->len + 3) / 4;
        size_t overhang = (numbytes * 4) - itr->len;
        assert(overhang < 4);

        size_t b = 0;
        char const *sb = tmpbuf;
        uint8_t *bytes = buf + byteoffset;
        uint8_t byte;

        while (b < numbytes)
        {
            byte = 0;
            int left = (b != numbytes-1? 4 : 4-overhang);

            for (int i = 0; i < left; ++i)
            {
                uint8_t code = DnaSeq::getcharcode(sb[i]);
                uint8_t shift = code << (6 - (2*i));
                byte |= shift;
            }

            bytes[b++] = byte;
            sb += 4;
        }

        offsets.push_back(byteoffset);
        byteoffset += b;
    }

    elapsed = MPI_Wtime() - elapsed;

    double mbspersecond = (totbases / 1048576.0) / elapsed;
    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << mbspersecond << " Mbs/second";
    logger.Flush("FASTA parsing rates (FastaIndex):");

    delete[] tmpbuf;
    delete[] readbuf;

    MPI_Exscan(&numreads, &firstid, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());
    if (commgrid->GetRank() == 0) firstid = 0;

    double avglen = static_cast<double>(totbases) / static_cast<double>(numreads);

    logger() << std::fixed << std::setprecision(2) << " my range [" << firstid << ".." << firstid+numreads << ") (" << numreads << " reads). ~" << avglen << " nts/read. (" << static_cast<double>(byteoffset) / (1024.0 * 1024.0) << " Mbs compressed)";
    logger.Flush("FASTA distributed among process ranks (FastaData):");
}

String FastaData::getsequence(size_t localid) const
{
    size_t len = readlens[localid];
    uint8_t const *bytes = buf + offsets[localid];
    Vector<char> s(len);

    for (size_t i = 0; i < len; ++i)
    {
        int code = (*bytes >> (6 - (2 * (i%4)))) & 3;
        s[i] = DnaSeq::getcodechar(code);

        if ((i+1) % 4 == 0)
            bytes++;
    }

    return String(s.begin(), s.end());
}

FastaData::FastaData(const FastaData& rhs) : commgrid(rhs.commgrid), offsets(rhs.offsets), readlens(rhs.readlens), firstid(rhs.firstid), bufsize(rhs.bufsize)
{
    buf = new uint8_t[rhs.bufsize];
    memcpy(buf, rhs.buf, bufsize);
}

FastaData FastaData::getrange(size_t pos, size_t count) const
{
    assert(pos + count <= numreads());

    FastaData rangefd(commgrid);
    rangefd.offsets.resize(count);
    rangefd.readlens.resize(count);

    auto lenitr = readlens.begin() + pos;
    std::copy(lenitr, lenitr + count, rangefd.readlens.begin());

    auto offsetitr = offsets.begin() + pos;
    size_t startoffset = *offsetitr;
    std::transform(offsetitr, offsetitr + count, rangefd.offsets.begin(), [&](size_t offset) { return offset - startoffset; });

    rangefd.bufsize = rangefd.offsets.back() + ((rangefd.readlens.back()+3)/4);
    rangefd.buf = new uint8_t[rangefd.bufsize];

    memcpy(rangefd.buf, rangefd.buf + offsets[pos], rangefd.bufsize);

    return rangefd;
}

FastaData& FastaData::operator+=(const FastaData& rhs)
{
    if (firstid + numreads() != rhs.firstid)
    {
        std::cerr << "Error: can only concatenate FastaData objects serving consecutive ranges" << std::endl;
        return *this;
    }

    size_t newbufsize = bufsize + rhs.bufsize;
    uint8_t *newbuf = new uint8_t[newbufsize];
    memcpy(newbuf, buf, bufsize);
    memcpy(newbuf + bufsize, rhs.buf, rhs.bufsize);

    offsets.reserve(numreads() + rhs.numreads());
    readlens.reserve(numreads() + rhs.numreads());

    for (size_t i = 0; i < rhs.numreads(); ++i)
    {
        readlens.push_back(rhs.readlens[i]);
        offsets.push_back(rhs.offsets[i] + bufsize);
    }

    delete[] buf;
    buf = newbuf;
    bufsize = newbufsize;

    return *this;
}

DistributedFastaData::DistributedFastaData(const FastaData& localdata) : commgrid(localdata.getcommgrid()), localdata(localdata)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    idoffsets.resize(nprocs);
    idoffsets[myrank] = localdata.getfirstid();

    MPI_Allgather(MPI_IN_PLACE, 1, MPI_SIZE_T, idoffsets.data(), 1, MPI_SIZE_T, comm);
}

FastaData DistributedFastaData::CollectRowReads() const
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    int myrowid = commgrid->GetRankInProcCol();
    int mycolid = commgrid->GetRankInProcRow();
    int numgridrows = commgrid->GetGridRows();

    size_t numreads = localdata.numreads();
    MPI_Allreduce(MPI_IN_PLACE, &numreads, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());

    assert(commgrid->GetGridRows() == commgrid->GetGridCols());

    size_t readspergridrow = numreads / numgridrows;
    size_t rowstartid = myrowid * readspergridrow;
    size_t colstartid = mycolid * readspergridrow;
    size_t rowendid = ((myrowid == numgridrows-1)? numreads : rowstartid + readspergridrow) - 1;
    size_t colendid = ((mycolid == numgridrows-1)? numreads : colstartid + readspergridrow) - 1;

    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << "row range ["    << rowstartid << ".." << rowendid << "] (" << (rowendid-rowstartid+1) << " reads). "
                                                   << "column range [" << colstartid << ".." << colendid << "] (" << (colendid-colstartid+1) << " reads). ";
    logger.Flush("Sequence grid distribution:");

    auto startitr = std::lower_bound(idoffsets.begin(), idoffsets.end(), rowstartid);
    auto enditr = std::lower_bound(startitr, idoffsets.end(), rowendid);


    return localdata;
}
