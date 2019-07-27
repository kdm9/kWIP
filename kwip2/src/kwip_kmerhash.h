#ifndef KWIP_KMERHASH_H_THACY8Z0
#define KWIP_KMERHASH_H_THACY8Z0

#include "kwip_config.h"

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
    uint64_t seed;
    bool canonicalise;
} kmer_iter_t;

void kmer_iter_init(kmer_iter_t *ctx, size_t k, uint64_t seed, bool canonicalise);
void kmer_iter_set_seq(kmer_iter_t *ctx, char *seq, size_t seqlen);

int kmer_iter_next_xxh(kmer_iter_t *ctx, uint64_t *hash);
int kmer_iter_next_nthash(kmer_iter_t *ctx, uint64_t *hash);

void kmer_iter_destroy(kmer_iter_t *ctx);

// One-shot kmer hash
uint64_t kmer_xxh(char *seq, size_t len, uint64_t seed, bool canonicalise);

#endif /* end of include guard: KWIP_KMERHASH_H_THACY8Z0 */
