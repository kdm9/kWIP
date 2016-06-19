#include "kmerhash.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

static const size_t N = 10000000;

#define XXHI

int main(int argc, char *argv[])
{
    //const char *seq = "ACGTACGTACGTNACGTNACGTACGT";
    const char *seq = "AAACAAGAATACCACGACTAGCAGGAGTATCATGATTCCCGCCTCGGCGTCTGCTTGGGTGTTTAA";
    const size_t k = 3;
    kmer_iter_t ctx;

    uint64_t hash = 0;

    const size_t c_sz = 1024;
    const size_t seed = 234;

    uint64_t *count = calloc(c_sz, sizeof(*count));

    for (size_t i = 0; i < N; i++) {
        kmer_iter_init(&ctx, seq, strlen(seq), k);
#ifdef XXHI
        while (kmer_iter_next_xxh(&ctx, &hash, seed)) {
#else
        while (kmer_iter_next(&ctx, &hash)) {
#endif
            count[hash % c_sz]++;
        }
    }

    printf("%llu\n", N);
    kmer_iter_init(&ctx, seq, strlen(seq), k);
#ifdef XXHI
    while (kmer_iter_next_xxh(&ctx, &hash, seed)) {
#else
    while (kmer_iter_next(&ctx, &hash)) {
#endif
        uint64_t cnt = count[hash % c_sz];
        if (cnt != N) {
            printf("h %llx c %llu\n", hash, cnt);
        }
    }

    free(count);
    return 0;
}
