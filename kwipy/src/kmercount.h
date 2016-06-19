#ifndef KMERCOUNT_H_VKHVXZLK
#define KMERCOUNT_H_VKHVXZLK
#include <stdlib.h>
#include <stdint.h>

#ifndef kc_eltype_t
    #define kc_eltype_t uint16_t
#endif

#ifndef KMERCOUNT_SAVE_CHUNKSIZE
    #define KMERCOUNT_SAVE_CHUNKSIZE 1048576
#endif
#define KMERCOUNT_USE_BLOSC
//#undef KMERCOUNT_USE_BLOSC

typedef struct {
    kc_eltype_t *cv;
    size_t len;
    size_t k;
    uint64_t seed;
} kmer_count_t;

void kmer_count_init(kmer_count_t *ctx, size_t len, size_t k, uint64_t seed);

kc_eltype_t kmer_count_count_h(kmer_count_t *ctx, uint64_t hash);
size_t kmer_count_count_s(kmer_count_t *ctx, const char *seq, size_t n);

kc_eltype_t kmer_count_get_h(kmer_count_t *ctx, uint64_t hash);

int kmer_count_save(kmer_count_t *ctx, const char *filename);
ssize_t kmer_count_consume_readfile(kmer_count_t *ctx, const char *filename);

void kmer_count_destroy(kmer_count_t *ctx);

#endif /* end of include guard: KMERCOUNT_H_VKHVXZLK */
