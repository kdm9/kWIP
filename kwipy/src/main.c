#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

#include <zlib.h>

#include "kmercount.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


int main(int argc, char *argv[])
{
    if (argc < 2) {
       fprintf(stderr, "USAGE: kmercount FASTX_FILE\n"); 
       return EXIT_FAILURE;
    }
    const size_t cvsize = 1000000000;

    gzFile fp = gzopen(argv[1], "r");
    kseq_t *seq = kseq_init(fp);
    kmer_count_t ctr;
    kmer_count_init(&ctr, cvsize, 21, 123);
    printf("Counting from %s\n", argv[1]);
    size_t total_kmers = 0;
    while(kseq_read(seq) >= 0) {
        total_kmers += kmer_count_count_s(&ctr, seq->seq.s, seq->seq.l);
    }
    printf("Done counting (%zu kmers)\n", total_kmers);

    kc_eltype_t max_c = 0;
    for (size_t i = 0; i < cvsize; i++) {
        kc_eltype_t this = ctr.cv[i];
        max_c = this > max_c ? this : max_c;
    }
    printf("max count is %u\n", (unsigned)max_c);
    kseq_destroy(seq);
    kmer_count_destroy(&ctr);
    gzclose(fp);
    return 0;
}
