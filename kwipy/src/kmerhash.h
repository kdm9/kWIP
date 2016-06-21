#ifndef KMERHASH_H_LDYCLTRK
#define KMERHASH_H_LDYCLTRK

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct {
    char *seq;
    char *rcseq;
    size_t len;
    size_t k;
    size_t i;
    uint64_t last_hash;
    uint64_t mask;
    bool canonicalise;
} kmer_iter_t;

void kmer_iter_init(kmer_iter_t *ctx, size_t k, bool canonicalise);
void kmer_iter_set_seq(kmer_iter_t *ctx, char *seq, size_t seqlen);

int kmer_iter_next(kmer_iter_t *ctx, uint64_t *hash);

int kmer_iter_next_xxh(kmer_iter_t *ctx, uint64_t *hash, uint64_t seed);

void kmer_iter_destroy(kmer_iter_t *ctx);

#endif /* end of include guard: KMERHASH_H_LDYCLTRK */
