#include "Buffer.h"
#include <assert.h>
#include <string.h>
#include <unistd.h>

buffer_t* init_buffer(size_t size)
{
    buffer_t *b = calloc(1, sizeof(*b));
    if (size == 0) return b;
    b->buf = calloc(size, 1);
    b->pos = b->len = 0;
    b->size = size;
    b->buf[b->len] = 0;
    return b;
}
size_t grow_buffer(buffer_t *b, size_t append_size)
{
    assert(b->size > 0);
    assert(append_size >= 0);

    size_t required_size = b->len + append_size + 1;

    if (required_size >= b->size)
    {
        while (required_size > b->size)
            b->size *= 2;

        assert(b->size >= required_size);

        b->buf = realloc(b->buf, b->size);
        assert(b->buf != NULL);
        memset(b->buf + b->len, 0, b->size - b->len); // initialize new memory to 0
    }

    assert(b->len + append_size < b->size);
    assert(b->buf);
    assert(b->size > 0);
    b->buf[b->size-1] = 0;
    return b->size - b->len;

}

char* get_start_buffer(buffer_t *b)
{
    return b->buf;
}

void free_buffer(buffer_t *b)
{
    if (!b) return;
    if (b->buf) free(b->buf);
    memset(b, 0, sizeof(buffer_t));
    free(b);
}
