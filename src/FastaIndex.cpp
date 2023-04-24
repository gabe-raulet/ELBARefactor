#include "FastaIndex.h"
#include "Logger.h"
#include <cstring>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>

typedef typename FastaIndex::faidx_record_t faidx_record_t;

faidx_record_t GetFaidxRecord(const String& line, Vector<String>& names)
{
    String name;
    faidx_record_t record;
    std::istringstream(line) >> name >> record.len >> record.pos >> record.bases;
    names.push_back(name);
    return record;
}

MPI_Count_type FastaIndex::get_idbalanced_partition(Vector<MPI_Count_type>& sendcounts)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    size_t numreads = GetNumRecords();

    sendcounts.resize(nprocs);

    MPI_Count_type readsperproc = numreads / nprocs;

    std::fill_n(sendcounts.begin(), nprocs-1, readsperproc);

    sendcounts.back() = numreads - (nprocs-1) * readsperproc;

    return numreads;
}

MPI_Count_type FastaIndex::get_membalanced_partition(Vector<MPI_Count_type>& sendcounts)
{
    return get_idbalanced_partition(sendcounts);
}

FastaIndex::FastaIndex(const String& fasta_fname, Grid commgrid, bool membalanced) : commgrid(commgrid), fasta_fname(fasta_fname)
{
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();

    Vector<MPI_Count_type> sendcounts; /* MPI_Scatterv sendcounts for faidx_record_t records (root only) */
    Vector<MPI_Displ_type> displs;     /* MPI_Scatterv displs for faidx_record_t records (root only)     */
    MPI_Count_type recvcount;          /* MPI_Scatterv recvcount for faidx_record_t records              */

    Vector<String> root_names;

    if (myrank == 0)
    {
        String line;
        std::ifstream filestream(GetFaidxFilename());

        while (std::getline(filestream, line))
        {
            records.push_back(GetFaidxRecord(line, root_names));
        }

        filestream.close();

        MPI_Count_type num_records;

        num_records = membalanced? get_membalanced_partition(sendcounts) : get_idbalanced_partition(sendcounts);

        displs.resize(nprocs);
        std::exclusive_scan(sendcounts.begin(), sendcounts.end(), displs.begin(), static_cast<MPI_Displ_type>(0));
    }

    /*
     * Root process tells each process how many faidx_record_t records it will be sent.
     */
    MPI_SCATTER(sendcounts.data(), 1, MPI_COUNT_TYPE, &recvcount, 1, MPI_COUNT_TYPE, 0, commgrid->GetWorld());

    myrecords.resize(recvcount);

    MPI_Datatype faidx_dtype_t;
    MPI_Type_contiguous(3, MPI_SIZE_T, &faidx_dtype_t);
    MPI_Type_commit(&faidx_dtype_t);
    MPI_SCATTERV(records.data(), sendcounts.data(), displs.data(), faidx_dtype_t, myrecords.data(), recvcount, faidx_dtype_t, 0, commgrid->GetWorld());
    MPI_Type_free(&faidx_dtype_t);
}

Vector<DnaSeq> FastaIndex::GetReadsFromRecords(const Vector<faidx_record_t>& records)
{
    Vector<DnaSeq> reads;

    size_t num_records = records.size();

    reads.reserve(num_records);

    const faidx_record_t& first_record = records.front();
    const faidx_record_t& last_record = records.back();

    MPI_Offset startpos = first_record.pos;
    MPI_Offset endpos = last_record.pos + last_record.len + (last_record.len / last_record.bases);

    MPI_File fh;
    MPI_File_open(commgrid->GetWorld(), fasta_fname.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    MPI_Offset filesize;
    MPI_File_get_size(fh, &filesize);

    endpos = std::min(endpos, filesize);

    MPI_Offset mychunksize = endpos - startpos;

    Vector<char> mychunk(mychunksize);

    MPI_FILE_READ_AT_ALL(fh, startpos, mychunk.data(), mychunksize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    size_t maxlen, totbases, offset = 0;

    MPI_Exscan(&num_records, &offset, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());

    maxlen   = std::accumulate(records.begin(), records.end(), static_cast<size_t>(0), [](size_t l, const faidx_record_t& rec) { return std::max(l, rec.len); });
    totbases = std::accumulate(records.begin(), records.end(), static_cast<size_t>(0), [](size_t l, const faidx_record_t& rec) { return l + rec.len; });

    char *seqbuf = new char[maxlen];

    double t0, t1;
    t0 = MPI_Wtime();
    for (auto itr = records.cbegin(); itr != records.cend(); ++itr)
    {
        size_t locpos = 0;
        ptrdiff_t chunkpos = itr->pos - startpos;
        ptrdiff_t remain = itr->len;
        char *bufptr = seqbuf;

        while (remain > 0)
        {
            size_t cnt = std::min(itr->bases, static_cast<size_t>(remain));
            std::memcpy(bufptr, &mychunk.data()[chunkpos + locpos], cnt);
            bufptr += cnt;
            remain -= cnt;
            locpos += (cnt+1);
        }

        reads.emplace_back(seqbuf, itr->len);
    }
    t1 = MPI_Wtime();

    double mbspersecond = (totbases / 1048576.0) / (t1-t0);
    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << mbspersecond << " Mbs/second";
    logger.Flush("FASTA parsing rates:");

    delete[] seqbuf;
    return reads;
}

Vector<DnaSeq> FastaIndex::GetMyReads()
{
    Vector<DnaSeq> myreads = GetReadsFromRecords(myrecords);

    size_t mynumreads = myreads.size();
    size_t mytotbases = std::accumulate(myrecords.begin(), myrecords.end(), static_cast<size_t>(0), [](size_t cur, const faidx_record_t& rec) { return cur + rec.len; });
    size_t mytotbytes = std::accumulate(myreads.begin(), myreads.end(), static_cast<size_t>(0), [](size_t cur, const auto& s) { return cur + s.numbytes(); });
    double myavglen = static_cast<double>(mytotbases) / static_cast<double>(mynumreads);

    size_t myreadoffset;
    MPI_Exscan(&mynumreads, &myreadoffset, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());
    if (commgrid->GetRank() == 0) myreadoffset = 0;

    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << " my range [" << myreadoffset << ".." << myreadoffset+mynumreads << "). ~" << myavglen << " nts/read. (" << static_cast<double>(mytotbytes) / (1024.0 * 1024.0) << " Mbs compressed)";
    logger.Flush("FASTA distributed among process ranks:");

    return myreads;
}
