#include "kwip_kmercount.h"

#define K 10
#define SKETCHSIZE 1048576
#define SEED 12345
#define CANONICAL true


void test_kmercount_counter(void **ctx)
{
    kmer_count_t ctr;
    char *seq = strdup("ACGTACGTAC");
    assert_non_null(seq);
    uint64_t khash = kmer_xxh(seq, strlen(seq), SEED, CANONICAL);
    assert_int_equal(khash, 0xcc5ba50198536bc8);

    kmer_count_init(&ctr, SKETCHSIZE, K, SEED, CANONICAL);

    kmer_count_count_h(&ctr, khash);
    assert_int_equal(kmer_count_get_h(&ctr, khash), 1);

    int ret = kmer_count_count_s(&ctr, seq, strlen(seq));
    assert_int_equal(ret, 0);
    assert_int_equal(kmer_count_get_h(&ctr, khash), 2);

    kmer_count_destroy(&ctr);
    free(seq);
}

void test_kmercount_loadsave(void **ctx)
{
    kmer_count_t ctr;
    int ret = 0;
    char *seq = strdup("ACGTACGTAC");
    uint64_t khash = kmer_xxh(seq, strlen(seq), SEED, CANONICAL);

    kmer_count_init(&ctr, SKETCHSIZE, K, SEED, CANONICAL);
    // Count kmer
    kmer_count_count_h(&ctr, khash);
    assert_int_equal(kmer_count_get_h(&ctr, khash), 1);

    // Save current state
    ret = kmer_count_save(&ctr, "counts.h5");
    assert_int_equal(ret, 0);
    // Count another to check that we can distinguish current state from loaded
    // state.
    kmer_count_count_h(&ctr, khash);
    assert_int_equal(kmer_count_get_h(&ctr, khash), 2);

    // Load counts
    ret = kmer_count_load(&ctr, "counts.h5");
    assert_int_equal(ret, 0);
    // Check state (should be 1, we saved before counting a 2nd time)
    assert_int_equal(kmer_count_get_h(&ctr, khash), 1);

    kmer_count_destroy(&ctr);
    remove("counts.h5");
    free(seq);
}

void test_kmercount_consume_file(void **ctx)
{
    const char *readfile = "data/10seq.fa";
    ssize_t ret = 0;
    kmer_count_t ctr;

    kmer_count_init(&ctr, SKETCHSIZE, K, SEED, CANONICAL);

    ret = kmer_count_consume_readfile(&ctr, readfile);
    assert_int_equal(ret, 0);
    assert_int_equal(ctr.num_reads, 10);
    assert_int_equal(ctr.num_kmers, 10 * (20 - K + 1));

    ret = kmer_count_consume_readfile(&ctr, "/no/file/exists/here");
    assert_int_equal(ret, 1);

    kmer_count_destroy(&ctr);
}


static const struct CMUnitTest suite_kmercount[] = {
    cmocka_unit_test(test_kmercount_counter),
    cmocka_unit_test(test_kmercount_consume_file),
    cmocka_unit_test(test_kmercount_loadsave),
};
