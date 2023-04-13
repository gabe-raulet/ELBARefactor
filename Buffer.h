#ifndef BUFFER_H_
#define BUFFER_H_

#include <stdlib.h>

#if defined (__cplusplus)
extern "C"
{
#endif

#ifndef MAX_ALLTOALL_MEM
#define MAX_ALLTOALL_MEM (128*1024*1024) /* 128 MB */
#endif

typedef struct
{
    char *buf;
    size_t pos, len, size;
} buffer_t;

buffer_t* init_buffer(size_t size);
size_t grow_buffer(buffer_t *buffer, size_t append_size);
char* get_start_buffer(buffer_t *buffer);
void free_buffer(buffer_t *buffer);

#endif
