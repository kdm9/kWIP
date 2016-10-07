#include <string.h>

#include "kwip_kmercount.h"

#define K 10
#define SKETCHSIZE 1048576
#define SEED 12345
#define CANONICAL true


TEST kwip_counter(void)
{
    kmer_count_t ctr;
    char *seq = strdup("ACGTACGTAC");
    uint64_t khash = kmer_xxh(seq, strlen(seq), SEED, CANONICAL);
    ASSERT_EQ(khash, 0xcc5ba50198536bc8);

    kmer_count_init(&ctr, SKETCHSIZE, K, SEED, CANONICAL);

    kmer_count_count_h(&ctr, khash);
    ASSERT_EQ(kmer_count_get_h(&ctr, khash), 1);

    size_t ret = kmer_count_count_s(&ctr, seq, strlen(seq));
    ASSERT_EQ(ret, 1);
    ASSERT_EQ(kmer_count_get_h(&ctr, khash), 2);

    kmer_count_destroy(&ctr);
    free(seq);
    PASS();
}

TEST kwip_loadsave(void)
{
    kmer_count_t ctr;
    int ret = 0;
    char *seq = strdup("ACGTACGTAC");
    uint64_t khash = kmer_xxh(seq, strlen(seq), SEED, CANONICAL);

    kmer_count_init(&ctr, SKETCHSIZE, K, SEED, CANONICAL);
    // Count kmer
    kmer_count_count_h(&ctr, khash);
    ASSERT_EQ(kmer_count_get_h(&ctr, khash), 1);

    // Save current state
    ret = kmer_count_save(&ctr, "counts.h5");
    ASSERT_EQ(ret, 0);
    // Count another to check that we can distinguish current state from loaded
    // state.
    kmer_count_count_h(&ctr, khash);
    ASSERT_EQ(kmer_count_get_h(&ctr, khash), 2);

    // Load counts
    ret = kmer_count_load(&ctr, "counts.h5");
    ASSERT_EQ(ret, 0);
    // Check state (should be 1, we saved before counting a 2nd time)
    ASSERT_EQ(kmer_count_get_h(&ctr, khash), 1);

    kmer_count_destroy(&ctr);
    remove("counts.h5");
    free(seq);
    PASS();
}

TEST consume_file(void)
{
    const char *readfile = "tests/data/10seq.fa";
    ssize_t ret = 0;
    kmer_count_t ctr;

    kmer_count_init(&ctr, SKETCHSIZE, K, SEED, CANONICAL);

    ret = kmer_count_consume_readfile(&ctr, readfile);
    ASSERT_EQ(ret, 10 /*reads*/ * (20 /*read length*/ - K + 1));

    ret = kmer_count_consume_readfile(&ctr, "/no/file/exists/here");
    ASSERT_EQ(ret, -1);

    kmer_count_destroy(&ctr);
    PASS();
}


SUITE(counting)
{
    RUN_TEST(counter);
    RUN_TEST(loadsave);
    RUN_TEST(consume_file);
}
