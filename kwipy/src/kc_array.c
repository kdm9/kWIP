#include <assert.h>

#include "kc_array.h"


#ifdef USE_BLOSC
#include "blosc/blosc_filter.h"
#endif

static const hsize_t chunkshape = CHUNKSIZE;
static const unsigned int compress_args[7] = {
    0, 0, 0, 0, /* 0 to 3 (inclusive) param slots are reserved. */
    9,          /* compression level */
    0,          /* 0: shuffle not active, 1: shuffle active */
    BLOSC_BLOSCLZ /* the compressor to use */
};

int
array_save(const char *filename, const char *dset_path, void *array,
           size_t items, hid_t dtype)
{
    int r;
    int return_code = 1;

    hid_t fid, sid, dset, plist = 0;

    assert(filename != NULL);
    assert(dset_path != NULL);
    assert(array != NULL);
    assert(items > 0);

    const hsize_t shape = items;

    // File
    fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(fid<0) goto end;

    // Dataspace
    sid = H5Screate_simple(1, &shape, &shape);
    if(sid<0) goto end;

    // Propertylist
    plist = H5Pcreate(H5P_DATASET_CREATE);
    if(plist<0) goto end;

    /* Chunked layout required for filters */
    r = H5Pset_chunk(plist, 1, &chunkshape);
    if(r<0) goto end;

#ifdef USE_BLOSC
    /* Register the filter with the library */
    r = register_blosc(NULL, NULL);
    if(r<0) goto end;

    /* Set the filter with 7 params */
    r = H5Pset_filter(plist, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, compress_args);
#else
    r = H5Pset_deflate(plist, 9);
#endif
    if(r<0) goto end;

    dset = H5Dcreate(fid, dset_path, dtype, sid, H5P_DEFAULT, plist,
                     H5P_DEFAULT);
    if(dset<0) goto end;

    r = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
    if(r<0) goto end;

    return_code = 0;

    end:
    if(dset>0)  H5Dclose(dset);
    if(sid>0)   H5Sclose(sid);
    if(plist>0) H5Pclose(plist);
    if(fid>0)   H5Fclose(fid);

    return return_code;
}


int
array_read(const char *filename, const char *dset_path, void **array, size_t *items, hid_t dtype)
{
    int r;
    int return_code = 1;

    hid_t fid, sid, dset, plist = 0;
    hsize_t shape;

    assert(filename != NULL);
    assert(dset_path != NULL);
    assert(array != NULL);
    assert(items != NULL);

#ifdef USE_BLOSC
    /* Register the filter with the library */
    r = register_blosc(NULL, NULL);
    if(r<0) goto end;
#endif

    // File
    fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(fid<0) goto end;

    // Dataset
    dset = H5Dopen(fid, dset_path, H5P_DEFAULT);
    if(dset<0) goto end;

    // Dataspace
    sid = H5Dget_space(dset);
    if(sid<0) goto end;

    // Shape
    int ndim = H5Sget_simple_extent_dims(sid, &shape, NULL);
    if (ndim != 1) goto end;
    if (shape > *items) {
        *array = realloc(*array, shape * H5Tget_size(dtype));
        assert(*array != NULL);
        *items = shape;
    }

    // Read array
    r = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, *array);
    if(r<0) goto end;

    return_code = 0;

end:
    if(dset>0)  H5Dclose(dset);
    if(sid>0)   H5Sclose(sid);
    if(plist>0) H5Pclose(plist);
    if(fid>0)   H5Fclose(fid);

    return return_code;
}
