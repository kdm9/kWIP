#ifndef KMERHASH_H_LDYCLTRK
#define KMERHASH_H_LDYCLTRK

#include <stdlib.h>
#include <stdint.h>

typedef struct {
    const char *seq;
    size_t len;
    size_t k;
    size_t i;
    uint64_t last_hash;
    uint64_t mask;
} kmer_iter_t;

void kmer_iter_init(kmer_iter_t *ctx, const char *seq, size_t n, size_t k);

int kmer_iter_next(kmer_iter_t *ctx, uint64_t *hash);

int kmer_iter_next_xxh(kmer_iter_t *ctx, uint64_t *hash, uint64_t seed);

void kmer_iter_destroy(kmer_iter_t *ctx);

#endif /* end of include guard: KMERHASH_H_LDYCLTRK */
