#include "SeqStore.h"

SeqStore::SeqStore(const FastaIndex& index) : commgrid(index.getcommgrid())
{
    auto myrecords = index.getmyrecords();
    size_t mynumrecords = myrecords.size();

    readlens.resize(mynumrecords);
    std::transform(myrecords.cbegin(), myrecords.cend(), readlens.begin(), [](const auto& record) { return record.len; });

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
    buf.resize(mybufsize);

    MPI_FILE_READ_AT_ALL(fh, startpos, buf.data(), mybufsize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}
