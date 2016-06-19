#include <assert.h>

#include "kmerhash.h"
#include "kmercount.h"
#ifdef KMERCOUNT_USE_BLOSC
#include "blosc/blosc_filter.h"
#endif

#include <hdf5.h>
#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

void
kmer_count_init(kmer_count_t *ctx, size_t len, size_t k, uint64_t seed)
{
    assert(ctx != NULL);
    assert(len > 0);
    assert(k > 1 && k <= 32);
    ctx->cv = calloc(len, sizeof(kc_eltype_t));
    assert(ctx->cv != NULL);
    ctx->len = len;
    ctx->k = k;
    ctx->seed = seed;
}

kc_eltype_t
kmer_count_count_h(kmer_count_t *ctx, uint64_t hash)
{
    assert(ctx != NULL);
    kc_eltype_t *bin = &ctx->cv[hash % ctx->len];
    kc_eltype_t cnt = *bin;
    return *bin = cnt + 1 < cnt ? cnt : cnt + 1;
}

size_t
kmer_count_count_s(kmer_count_t *ctx, const char *seq, size_t n)
{
    kmer_iter_t itr;
    kmer_iter_init(&itr, seq, n, ctx->k);
    uint64_t hash = 0;
    size_t num_kmers = 0;
    for (; kmer_iter_next_xxh(&itr, &hash, ctx->seed); num_kmers++) {
        kmer_count_count_h(ctx, hash);
    }
    return num_kmers;
}

kc_eltype_t
kmer_count_get_h(kmer_count_t *ctx, uint64_t hash)
{
    assert(ctx != NULL);
    return ctx->cv[hash % ctx->len];
}

int
kmer_count_save(kmer_count_t *ctx, const char *filename)
{
    const hsize_t chunkshape = KMERCOUNT_SAVE_CHUNKSIZE;
    char *version, *date;
    unsigned int compress_args[7];
    int return_code = 1;

    hid_t fid, sid, dset, plist = 0;

    assert(ctx != NULL);
    const hsize_t shape = ctx->len;

    int r;

    // File
    fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(fid<0) goto failed;

    // dataspace
    sid = H5Screate_simple(1, &shape, &shape);
    if(sid<0) goto failed;


    // Propertylist
    plist = H5Pcreate(H5P_DATASET_CREATE);
    if(plist<0) goto failed;

    /* Chunked layout required for filters */
    r = H5Pset_chunk(plist, 1, &chunkshape);
    if(r<0) goto failed;

#ifdef KMERCOUNT_USE_BLOSC
    /* Register the filter with the library */
    r = register_blosc(&version, &date);
    if(r<0) goto failed;

    /* 0 to 3 (inclusive) param slots are reserved. */
    compress_args[4] = 9;       /* compression level */
    compress_args[5] = 0;       /* 0: shuffle not active, 1: shuffle active */
    compress_args[6] = BLOSC_BLOSCLZ; /* the actual compressor to use */

    /* Set the filter with 7 params */
    r = H5Pset_filter(plist, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, compress_args);
#else
    r = H5Pset_deflate(plist, 9);
#endif
    if(r<0) goto failed;

    dset = H5Dcreate(fid, "count", H5T_STD_U16LE, sid, H5P_DEFAULT, plist,
                     H5P_DEFAULT);

    if(dset<0) goto failed;

    r = H5Dwrite(dset, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, ctx->cv);
    if(r<0) goto failed;

    /*
    r = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data_out);
    if(r<0) goto failed;

    for(i=0;i<SIZE;i++){
        if(data[i] != data_out[i]) goto failed;
    }
    */

    return_code = 0;

    failed:

    if(dset>0)  H5Dclose(dset);
    if(sid>0)   H5Sclose(sid);
    if(plist>0) H5Pclose(plist);
    if(fid>0)   H5Fclose(fid);

    return return_code;
}

ssize_t
kmer_count_consume_readfile(kmer_count_t *ctx, const char *filename)
{
    gzFile fp = gzopen(filename, "r");
    if (fp == NULL) return -1;
    kseq_t *seq = kseq_init(fp);
    if (seq == NULL) return -1;

    size_t num_reads = 0;
    while(kseq_read(seq) >= 0) {
        kmer_count_count_s(ctx, seq->seq.s, seq->seq.l);
        num_reads++;
    }
    kseq_destroy(seq);
    gzclose(fp);
    return num_reads;
}

void
kmer_count_destroy(kmer_count_t *ctx)
{
    if (ctx != NULL) {
        if (ctx->cv != NULL) {
            free(ctx->cv);
        }
    }
}
