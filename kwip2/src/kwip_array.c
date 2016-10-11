#include "kwip_array.h"

#include <assert.h>



#ifdef USE_BLOSC
#include "blosc_filter.h"
static const unsigned int compress_args[7] = {
    0, 0, 0, 0, /* 0 to 3 (inclusive) param slots are reserved. */
    9,          /* compression level */
    1,          /* 0: shuffle not active, 1: shuffle active */
    BLOSC_BLOSCLZ /* the compressor to use */
};
#endif

int
array_save(const char *filename, const char *dset_path, void *array,
           size_t items, hid_t dtype)
{
    int r;
    int return_code = 1;
    const hsize_t chunkshape = KWIP_CHUNKSIZE;

    hid_t fid = 0, sid = 0, dset = 0, plist = 0;

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
    r = H5Pset_deflate(plist, 1);
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

    hid_t fid = 0, sid = 0, dset = 0, plist = 0;
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


int
array_blockiter_init(array_blockiter_t *ctx, const char *filename,
                     const char *dset_path)
{
    assert(filename != NULL);
    assert(dset_path != NULL);
    ctx->chunknum = 0;

    int r;
    int return_code = 1;
#ifdef USE_BLOSC
    /* Register the filter with the library */
    r = register_blosc(NULL, NULL);
    if(r<0) goto end;
#endif

    ctx->h5f = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (ctx->h5f < 0) goto end;

    // Dataset
    ctx->dset = H5Dopen(ctx->h5f, dset_path, H5P_DEFAULT);
    if (ctx->dset < 0) goto end;

    // Dataspace
    ctx->dspace = H5Dget_space(ctx->dset);
    if (ctx->dspace < 0) goto end;

    // Datatype
    ctx->dtype = H5Dget_type(ctx->dset);
    if (ctx->dspace < 0) goto end;

    // Shape
    int ndim = H5Sget_simple_extent_dims(ctx->dspace, &ctx->shape, NULL);
    if (ndim != 1) goto end;

    ctx->num_chunks = (ctx->shape + (KWIP_CHUNKSIZE / 2)) / KWIP_CHUNKSIZE;

    return_code = 0;
end:
    return return_code;
}

inline int array_blockiter_done(array_blockiter_t *ctx)
{
    if (ctx == NULL) return -1;
    if (ctx->chunknum >= ctx->num_chunks) {
        return 1;
    }
    return 0;
}

int
array_blockiter_next(array_blockiter_t *ctx, void **block, size_t *blocklen,
                     size_t maxsz)
{
    herr_t res;
    assert(ctx != NULL);
    assert(block != NULL);
    assert(blocklen != NULL);
    assert(maxsz >= KWIP_CHUNKSIZE);

    if (array_blockiter_done(ctx)) {
        // Someone called next on an iterator that has already finished
        return 0;
    }

    if (*block == NULL || *blocklen < KWIP_CHUNKSIZE) {
        *block = realloc(*block, KWIP_CHUNKSIZE * H5Tget_size(ctx->dtype));
        if (*block == NULL) {
            fprintf(stderr, "realloc failed (OUT OF MEMORY?)\n");
            return -1;
        }
    }


    // setup dataspace & hyperslab for read
    hid_t space = H5Dget_space(ctx->dset);
    assert(space != 0);
    hsize_t start = ctx->chunknum * KWIP_CHUNKSIZE;
    hsize_t to_end = ctx->shape - start;
    hsize_t count = to_end > KWIP_CHUNKSIZE ? KWIP_CHUNKSIZE : to_end;
    hsize_t capacity = KWIP_CHUNKSIZE;
    hid_t memspace = H5Screate_simple(1, &capacity, NULL);
    assert(memspace != 0);
    
    res = H5Sselect_hyperslab(space, H5S_SELECT_SET, &start, NULL, &count, NULL);
    if (res < 0) {
        fprintf(stderr, "Error selecting hyperslab\n");
        return -1;
    }

    res = H5Dread(ctx->dset, ctx->dtype, memspace, space, H5P_DEFAULT, *block);
    if (res < 0) {
        fprintf(stderr, "Error reading chunk\n");
        return -1;
    }
    *blocklen = count;
    ctx->chunknum++;
    return 0;
}

void
array_blockiter_destroy(array_blockiter_t *ctx)
{
    if (ctx->dtype > 0) H5Tclose(ctx->dtype);
    if (ctx->dspace > 0) H5Sclose(ctx->dspace);
    if (ctx->dset > 0) H5Dclose(ctx->dset);
    if (ctx->h5f > 0) H5Fclose(ctx->h5f);
    ctx->chunknum = 0;
    ctx->shape = 0;
}
