#ifndef ARRAYIO_H_6AYJOGQ7
#define ARRAYIO_H_6AYJOGQ7

#include <hdf5.h>
#ifndef CHUNKSIZE
    #define CHUNKSIZE 1048576
#endif
#define USE_BLOSC


int array_save(const char *filename, const char *dset_path, void *array,
               size_t items, hid_t dtype);
int array_read(const char *filename, const char *dset_path, void **array,
               size_t *items, hid_t dtype);


#endif /* end of include guard: ARRAYIO_H_6AYJOGQ7 */
