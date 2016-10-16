#ifndef KWIP_ARRAY_H_J1IUDTGQ
#define KWIP_ARRAY_H_J1IUDTGQ

#include "kwip_config.h"

#include <hdf5.h>


int array_save(const char *filename, const char *dset_path, void *array,
               size_t items, hid_t dtype);
int array_read(const char *filename, const char *dset_path, void **array,
               size_t *items, hid_t dtype);

typedef struct {
    hid_t h5f;
    hid_t dspace;
    hid_t dtype;
    hid_t dset;
    size_t num_chunks;
    hsize_t chunknum;
    hsize_t shape;
} array_blockiter_t;

int array_blockiter_init(array_blockiter_t *ctx, const char *filename,
                         const char *dset_path);
// array_blockiter_done returns -1 on error, 1 if done, 0 otherwise
extern int array_blockiter_done(array_blockiter_t *ctx);
int array_blockiter_next(array_blockiter_t *ctx, void **block,
                         size_t *blocklen, size_t maxsize);
void array_blockiter_destroy(array_blockiter_t *ctx);

#endif /* end of include guard: KWIP_ARRAY_H_J1IUDTGQ */
