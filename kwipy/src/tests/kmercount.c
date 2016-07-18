#include <string.h>
#include "kmercount.h"

#define K 10
#define NK 1048576
#define SEED 12345
#define CANONICAL true



TEST counter(void)
{
    kmer_count_t ctr;
    size_t ret = 0;
    char *seq = strdup("ACGTACGTAC");
    uint64_t khash = kmer_xxh(seq, strlen(seq), SEED, CANONICAL);
    ASSERT_EQ(khash, 0xcc5ba50198536bc8);

    kmer_count_init(&ctr, NK, K, SEED, CANONICAL);

    kmer_count_count_h(&ctr, khash);
    ASSERT_EQ(kmer_count_get_h(&ctr, khash), 1);

    ret = kmer_count_count_s(&ctr, seq, strlen(seq));
    ASSERT_EQ(ret, 1);
    ASSERT_EQ(kmer_count_get_h(&ctr, khash), 2);

    kmer_count_destroy(&ctr);
    free(seq);
    PASS();
}


SUITE(counting)
{
    RUN_TEST(counter);
}
