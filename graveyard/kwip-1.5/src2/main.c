#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>

#include <zlib.h>

#include "kmercount.h"


int main(int argc, char *argv[])
{
    if (argc < 3) {
       fprintf(stderr, "USAGE: kmercount COUNTS_FILE FASTX_FILE\n");
       return EXIT_FAILURE;
    }
    const size_t cvsize = 1000000000;
    const char *savefile = argv[1];

    clock_t start;
    double secs;
    kmer_count_t ctr;
    kmer_count_init(&ctr, cvsize, 20, 123, true);
    for (int fidx = 2; fidx < argc; fidx++) {
        const char *readfile = argv[fidx];
        start = clock();
        printf("Counting from '%s'...\n", readfile);
        size_t nreads = kmer_count_consume_readfile(&ctr, readfile);
        secs = (double)(clock() - start) / CLOCKS_PER_SEC;
        printf("\t- done! (%zu reads, %0.2fs, %0.1fK r/s)\n", nreads,
                secs, (double)(nreads / 1000) / secs);
    }

    kc_eltype_t max_c = 0;
    for (size_t i = 0; i < cvsize; i++) {
        kc_eltype_t this = ctr.cv[i];
        max_c = this > max_c ? this : max_c;
    }
    printf("max count is %u\n", (unsigned)max_c);
    start = clock();
    kmer_count_save(&ctr, savefile);
    secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("Wrote '%s' in %0.1f sec\n", savefile, secs);
    kmer_count_destroy(&ctr);


    kmer_count_load(&ctr, savefile);
    max_c = 0;
    for (size_t i = 0; i < cvsize; i++) {
        kc_eltype_t this = ctr.cv[i];
        max_c = this > max_c ? this : max_c;
    }
    printf("max count is (still) %u\n", (unsigned)max_c);
    kmer_count_destroy(&ctr);
    return 0;
}
