#include "kwip_config.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>

#include <zlib.h>
#include <getopt.h>

#include <clogged/clogged.h>

#include "kwip_utils.h"

#include "kwip_kmercount.h"


/*******************************************************************************
*                                    COUNT                                    *
*******************************************************************************/

static void count_help(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    kwip count [options] COUNT_FILE SEQUENCE_FILE ...\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -k, --ksize        K-mer length [default: 20]\n");
    fprintf(stream, "    -s, --sketchsize   Sketch size, in bins [default: 100M]\n");
    //fprintf(stream, "    -v, --verbose      Be verbose (default)\n");
    //fprintf(stream, "    -V, --version      Print version\n");
    //fprintf(stream, "    -q, --quiet        Be quiet\n");
}

static const char count_optstring[] = "k:s:vVqh";
static const struct option count_options[] = {
    {"ksize"      , required_argument , NULL , 'k'} ,
    {"sketchsize" , required_argument , NULL , 's'} ,
    {"verbose"    , no_argument       , NULL , 'v'} ,
    {"version"    , no_argument       , NULL , 'V'} ,
    {"quiet"      , no_argument       , NULL , 'q'} ,
    {"help"       , no_argument       , NULL , 'h'} ,
    {NULL         , 0                 , NULL , 0}
};

int count_main(int argc, char *argv[])
{
    size_t sketchsize = 100<<20;
    size_t ksize = 0;
    clg_logger_t *log = clg_logger_create();
    clg_logger_default(log, CLG_LOG_PROGRESS);

    int c = 0;
    while (1) {
        int option_index = 0;
        c = getopt_long(argc, argv, count_optstring, count_options, &option_index);
        if (c < 0) break;

        switch(c) {
            case 's':
                sketchsize = kwip_parse_size(optarg);
                break;
            case 'k':
                ksize = strtol(optarg, NULL, 10);
                break;
            case 'V':
                break;
            case 'v':
                clg_logger_set_level(log, CLG_LOG_DEBUG);
                break;
            case 'q':
                clg_logger_set_level(log, CLG_LOG_INFO);
                break;
            case 'h':
                count_help(stdout);
                return EXIT_SUCCESS;
            case '?':
                /* getopt_long will have already printed an error */
                break;
            default:
                count_help(stderr);
                return EXIT_FAILURE;
        }
    }

    if (sketchsize < 1<<10) {
        fprintf(stderr, "ERROR: Sketch size far too small. Increase to at least 1M.\n\n");
        count_help(stderr);
        return EXIT_FAILURE;
    }
    if (ksize < 10 || ksize > 32) {
        fprintf(stderr, "ERROR: kmer size must be 10 <= k <= 32.\n\n");
        count_help(stderr);
        return EXIT_FAILURE;
    }
    if (optind >= argc) {
        fprintf(stderr, "ERROR: must provide Output and Input file names\n\n");
        count_help(stderr);
        return EXIT_FAILURE;
    }

    const char *savefile = argv[optind++];

    kmer_count_t ctr;
    uint64_t seed = 1234;
    bool canonicalise = true;
    kmer_count_init(&ctr, sketchsize, ksize, seed, canonicalise);
    kmer_count_set_logger(&ctr, log);

    clock_t start;
    double secs;
    for (; optind < argc; optind++) {
        const char *readfile = argv[optind];
        start = clock();
        clg_log_fmt_info(log, "Counting from '%s' ...\n", readfile);
        size_t nreads = kmer_count_consume_readfile(&ctr, readfile);
        secs = (double)(clock() - start) / CLOCKS_PER_SEC;
        clg_log_fmt_info(log, "\t- done! (%zu reads, %0.2fs, %0.1fK r/s)\n", nreads,
                secs, (double)(nreads / 1000) / secs);
    }

    kmer_count_save(&ctr, savefile);
    secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    clg_log_fmt_info(log, "Wrote '%s' in %0.1f sec\n", savefile, secs);
    kmer_count_destroy(&ctr);
    clg_logger_destroy(log);
    optind = 1;
    return 0;
}

void global_help(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    kwip SUBCOMMAND [options]\n");
    fprintf(stream, "\n");
    fprintf(stream, "Where subcommand is one of:\n");
    fprintf(stream, "\t- count: Count kmers into sketches\n");
}

int main(int argc, char *argv[])
{
    if (argc < 2 || \
            strcmp(argv[1], "-h") == 0 || \
            strcmp(argv[1], "--help") == 0) {
        global_help(stdout);
        return EXIT_SUCCESS;
    }

    if (strcmp(argv[1], "count") == 0) {
        return count_main(argc - 1, argv + 1);
    } else {
        global_help(stderr);
        return EXIT_FAILURE;
    }
}
