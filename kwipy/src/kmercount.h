#ifndef KMERCOUNT_H_VKHVXZLK
#define KMERCOUNT_H_VKHVXZLK


#define kc_eltype_t uint16_t

typedef struct {
    kc_eltype_t *cv;
    size_t len;
    size_t k;
} kmer_count_t;

void kmer_count_init(kmer_count_t *ctx, size_t len, size_t k);

kc_eltype_t kmer_count_count_h(kmer_count_t *ctx, uint64_t hash);
void kmer_count_count_s(kmer_count_t *ctx, const char *seq, size_t n, uint64_t seed);

kc_eltype_t kmer_count_get_h(kmer_count_t *ctx, uint64_t hash);

void kmer_count_destroy(kmer_count_t *ctx);

#endif /* end of include guard: KMERCOUNT_H_VKHVXZLK */
