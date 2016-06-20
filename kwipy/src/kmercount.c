#include <assert.h>

#include "kmerhash.h"
#include "kmercount.h"
#include "kc_array.h"

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
    return array_save(filename, "counts", ctx->cv, ctx->len, H5T_NATIVE_UINT16);
}

int
kmer_count_load(kmer_count_t *ctx, const char *filename)
{
    return array_read(filename, "counts", (void **)&ctx->cv, &ctx->len, H5T_NATIVE_UINT16);
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
            ctx->len = 0;
            ctx->cv = NULL;
        }
    }
}
