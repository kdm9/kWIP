#include "kwip_kmercount.h"

#include <assert.h>

#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

void
kmer_count_init(kmer_count_t *ctx, size_t cvlen, size_t k, uint64_t seed, bool canonicalise)
{
    assert(ctx != NULL);
    assert(cvlen > 0);
    assert(k >= 0 && k <= 32);
    ctx->cv = calloc(cvlen, sizeof(kc_eltype_t));
    assert(ctx->cv != NULL);
    ctx->cvlen = cvlen;
    ctx->k = k;
    kmer_iter_init(&ctx->itr, ctx->k, seed, canonicalise);
}

void
kmer_count_set_logger(kmer_count_t *ctx, clg_logger_t *log)
{
    assert(ctx != NULL);
    ctx->log = log;
}

kc_eltype_t
kmer_count_count_h(kmer_count_t *ctx, uint64_t hash)
{
    assert(ctx != NULL);
    size_t idx = hash % ctx->cvlen;
    unsigned cnt = ctx->cv[idx];
    cnt = cnt == ((kc_eltype_t) -1) ? cnt : cnt + 1;
    ctx->cv[idx] = cnt;
    return cnt;
}

size_t
kmer_count_count_s(kmer_count_t *ctx, char *seq, size_t n)
{
    size_t num_kmers = 0;
    kmer_iter_set_seq(&ctx->itr, seq, n);
    for (uint64_t hash = 0; kmer_iter_next_xxh(&ctx->itr, &hash); num_kmers++) {
        kmer_count_count_h(ctx, hash);
    }
    return num_kmers;
}

kc_eltype_t
kmer_count_get_h(kmer_count_t *ctx, uint64_t hash)
{
    assert(ctx != NULL);
    return ctx->cv[hash % ctx->cvlen];
}

int
kmer_count_save(kmer_count_t *ctx, const char *filename)
{
    return array_save(filename, "counts", ctx->cv, ctx->cvlen, H5T_KMERCOUNT);
}

int
kmer_count_load(kmer_count_t *ctx, const char *filename)
{
    return array_read(filename, "counts", (void **)&ctx->cv, &ctx->cvlen, H5T_KMERCOUNT);
}

static inline ssize_t
kmer_count_consume_fp(kmer_count_t *ctx, gzFile fp)
{
    if (fp == NULL) return -1;
    kseq_t *seq = kseq_init(fp);
    if (seq == NULL) return -1;

    size_t num_reads = 0;
    size_t num_kmers = 0;
    while(kseq_read(seq) >= 0) {
        num_kmers += kmer_count_count_s(ctx, seq->seq.s, seq->seq.l);
        num_reads++;
        if (ctx->log != NULL && num_reads % 200000 == 0) {
            clg_log_fmt_progress(ctx->log, "\t- Counted %0.1fM reads ...\n", (float)num_reads / 1000000);
        }
    }
    kseq_destroy(seq);
    return num_reads;
}

ssize_t
kmer_count_consume_readfile(kmer_count_t *ctx, const char *filename)
{
    gzFile fp = gzopen(filename, "r");
    ssize_t num_reads = kmer_count_consume_fp(ctx, fp);
    gzclose(fp);
    return num_reads;
}

ssize_t
kmer_count_consume_fd(kmer_count_t *ctx, int fd)
{
    gzFile fp = gzdopen(fd, "r");
    ssize_t num_reads = kmer_count_consume_fp(ctx, fp);
    gzclose(fp);
    return num_reads;
}


void
kmer_count_destroy(kmer_count_t *ctx)
{
    if (ctx != NULL) {
        if (ctx->cv != NULL) {
            free(ctx->cv);
            ctx->cvlen = 0;
            ctx->cv = NULL;
        }
        kmer_iter_destroy(&ctx->itr);
    }
}
