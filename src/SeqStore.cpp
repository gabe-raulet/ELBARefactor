#include "SeqStore.h"
#include "Logger.h"
#include <cstring>
#include <cassert>

SeqStore::SeqStore(const FastaIndex& index) : commgrid(index.getcommgrid())
{
    auto myrecords = index.getmyrecords();
    size_t mynumrecords = myrecords.size();

    offsets.reserve(mynumrecords);
    readlens.reserve(mynumrecords);

    const auto& first_record = myrecords.front();
    const auto& last_record = myrecords.back();

    MPI_Offset startpos = first_record.pos;
    MPI_Offset endpos = last_record.pos + last_record.len + (last_record.len / last_record.bases);

    MPI_File fh;
    MPI_File_open(commgrid->GetWorld(), index.GetFastaFilename().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    MPI_Offset filesize;
    MPI_File_get_size(fh, &filesize);

    if (endpos > filesize) endpos = filesize;

    MPI_Offset mybufsize = endpos - startpos;

    buf = new char[mybufsize];

    MPI_FILE_READ_AT_ALL(fh, startpos, buf, mybufsize, MPI_CHAR, MPI_STATUS_IGNORE);

    MPI_File_close(&fh); /* implicit barrier */
    double elapsed = MPI_Wtime();

    char *rptr;
    char *wptr;

    for (auto itr = myrecords.cbegin(); itr != myrecords.cend(); ++itr)
    {
        readlens.push_back(itr->len);
        offsets.push_back(itr->pos - startpos); /* itr->pos is global file offset. subtract local process file start position to get my offset */

        wptr = rptr = buf + offsets.back();

        for (size_t i = 0; i < itr->len; ++i)
        {
            *wptr++ = DnaSeq::getcharchar(*rptr++);

            if ((i+1) % itr->bases == 0) /* '\n' character ever itr->bases. increment read pointer */
                rptr++;
        }

        *wptr = '\0';

        assert(strlen(buf + offsets.back()) == itr->len);
    }

    elapsed = MPI_Wtime() - elapsed;

    double worktime, spantime, costtime;
    MPI_Reduce(&elapsed, &worktime, 1, MPI_DOUBLE, MPI_SUM, 0, commgrid->GetWorld());
    MPI_Reduce(&elapsed, &spantime, 1, MPI_DOUBLE, MPI_MAX, 0, commgrid->GetWorld());

    costtime = spantime * commgrid->GetSize();

    if (commgrid->GetRank() == 0)
    {
        std::cout << std::fixed << std::setprecision(3) << "[" << worktime << ", " << costtime << ", " << spantime << "] [work, cost, span] (processor seconds)" << std::endl;
    }

    MPI_Barrier(commgrid->GetWorld());
}

SeqStore::~SeqStore()
{
    delete[] buf;
}

size_t SeqStore::GetNumReads() const
{
    size_t mynumreads = GetMyNumReads();
    MPI_Allreduce(MPI_IN_PLACE, &mynumreads, 1, MPI_SIZE_T, MPI_SUM, commgrid->GetWorld());
    return mynumreads;
}

Vector<size_t> SeqStore::GetGlobalReadCounts() const
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    Vector<size_t> counts(nprocs);
    counts[myrank] = GetMyNumReads();

    MPI_Allgather(MPI_IN_PLACE, 1, MPI_SIZE_T, counts.data(), 1, MPI_SIZE_T, comm);

    return counts;
}

Vector<size_t> SeqStore::GetGlobalReadOffsets() const
{
    Vector<size_t> counts = GetGlobalReadCounts();
    Vector<size_t> readoffsets(counts.size());
    std::exclusive_scan(counts.begin(), counts.end(), readoffsets.begin(), static_cast<size_t>(0));
    return readoffsets;
}
