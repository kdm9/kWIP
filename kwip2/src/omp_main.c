#include "kwip_config.h"

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>

#include <zlib.h>
#include <getopt.h>

#include <clogged/clogged.h>

#include "kwip_utils.h"

#include "kwip_kmercount.h"
#include "kwip_kernelcalc.h"
#include "kwip_omp_kernel.h"
#include "kwip_metric_ip.h"


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
    fprintf(stream, "    -v, --verbose      Be verbose (default)\n");
    fprintf(stream, "    -V, --version      Print version\n");
    fprintf(stream, "    -q, --quiet        Be quiet\n");
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
    int res = 0;
    for (; optind < argc; optind++) {
        const char *readfile = argv[optind];
        start = clock();
        clg_log_fmt_info(log, "Counting from '%s' ...\n", readfile);
        res = kmer_count_consume_readfile(&ctr, readfile);
        if (res != 0) {
            // TODO: error
        }
        secs = (double)(clock() - start) / CLOCKS_PER_SEC;
        float readpersec = (double)(ctr.num_reads / 1000) / secs;
        clg_log_fmt_info(log, "\t- done! (%zu reads, %0.2fs, %0.1fK r/s)\n",
                         ctr.num_reads, secs, readpersec);
    }

    start = clock();
    res = kmer_count_save(&ctr, savefile);
    if (res != 0) {
        // TODO: error
    }
    secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    clg_log_fmt_info(log, "Wrote '%s' in %0.1f sec\n", savefile, secs);
    kmer_count_destroy(&ctr);
    clg_logger_destroy(log);
    optind = 1;
    return 0;
}


/*******************************************************************************
*                                    kernel                                    *
*******************************************************************************/

static void kernel_help(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    kwip kernel [options] COUNT_FILE ...\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -m, --metric   Metric for kernel calculation. Must be one of\n");
    fprintf(stream, "                   wip, ip or manhattan. [default: ip]\n"); //TODO: change to wip
    fprintf(stream, "    -d, --outdir   Output directory (REQUIRED) \n");
    fprintf(stream, "    -v, --verbose  Be verbose (default)\n");
    fprintf(stream, "    -q, --quiet    Be quiet\n");
}

static const char kernel_optstring[] = "m:o:t:vqh";
static const struct option kernel_options[] = {
    {"threads" , required_argument , NULL , 't'},
    {"metric"  , required_argument , NULL , 'm'},
    {"outdir"  , required_argument , NULL , 'o'},
    {"verbose" , no_argument       , NULL , 'v'},
    {"quiet"   , no_argument       , NULL , 'q'},
    {"help"    , no_argument       , NULL , 'h'},
    {NULL      , 0                 , NULL , 0}
};

int kernel_main(int argc, char *argv[])
{
    char *metric = strdup("ip");
    char *outdir = strdup("kernel_out");
    int threads = 1;
    int retval = 0;
    clg_logger_t *log = clg_logger_create();
    clg_logger_default(log, CLG_LOG_PROGRESS);

    int c = 0;
    while (1) {
        int option_index = 0;
        c = getopt_long(argc, argv, kernel_optstring, kernel_options, &option_index);
        if (c < 0) break;

        switch(c) {
            case 't':
                threads = atoi(optarg);
                break;
            case 'm':
                metric = strdup(optarg);
                break;
            case 'o':
                outdir = strdup(optarg);
                break;
            case 'v':
                clg_logger_set_level(log, CLG_LOG_DEBUG);
                break;
            case 'q':
                clg_logger_set_level(log, CLG_LOG_INFO);
                break;
            case 'h':
                kernel_help(stdout);
                retval = 0; goto exit;
            case '?':
                break;
            default:
                kernel_help(stderr);
                retval = -1; goto exit;
        }
    }

    int res = 0;

    res = kwip_mkdirp(outdir);
    if (res != 0) {
        clg_log_fmt_error(log, "Failed to make output directory '%s': %s\n",
                          outdir, strerror(errno));
        retval = -1; goto exit;
    }

    if (optind + 1 >= argc) {
        clg_log_msg_error(log,
                "ERROR: must provide at least two sketch files to compare\n\n");

        kernel_help(stderr);
        retval = -1; goto exit;
    }

    kwip_kerncalc_t ctx;
    res = kerncalc_init(&ctx);
    if (res != 0) {
        clg_log_msg_error(log, "ERROR: Could not initialise kernel calculator.\n");
        retval = -1; goto exit;
    }

    if (strcasecmp(metric, "ip") == 0) {
        // TODO add normalisation etc. here if we use that
        res = kerncalc_set_kernelfunction(&ctx, metric_ip_kernel);
        if (res != 0) {
            clg_log_msg_error(log, "ERROR: Can't set kernel function.\n");
            retval = -1; goto exit;
        }
    } else {
        if (res != 0) {
            clg_log_fmt_error(log, "ERROR: Unknown metric '%s'\n", metric);
            retval = -1; goto exit;
        }
    }

    res = kerncalc_set_checkpoint_dir(&ctx, outdir);
    if (res != 0) {
        clg_log_msg_error(log, "ERROR: Can't set output directory.\n");
        retval = -1; goto exit;
    }



    for (; optind < argc; optind++) {
        res = kerncalc_add_sample(&ctx, argv[optind], NULL);
        if (res != 0) {
            clg_log_fmt_error(log, "ERROR: Failed to add sample '%s'.\n", argv[optind]);
            retval = -1; goto exit;
        }
    }
    res = kerncalc_finalise(&ctx, NULL, NULL);
    if (res != 0) {
        clg_log_msg_error(log, "ERROR: Failed to finalise kernel calculator.\n");
        retval = -1; goto exit;
    }

    clg_log_fmt_info(log, "Starting pairwise similiarity calculation using %d threads\n",
                     threads);
    res = kerncalc_pairwise_omp(&ctx, log, threads);
    if (res != 0) {
        clg_log_msg_error(log, "ERROR: Pairwise comparison failed.\n");
        retval = -1; goto exit;
    }
    clg_log_msg_info(log, "Completed pairwise similiarity calculation\n");

    res = kerncalc_save(&ctx);
    if (res != 0) {
        clg_log_msg_error(log, "ERROR: saving kernel values failed\n");
        retval = -1; goto exit;
    }

exit:
    kerncalc_destroy(&ctx);
    clg_logger_destroy(log);
    kwip_free(outdir);
    kwip_free(metric);
    optind = 1;
    return retval;
}

/*******************************************************************************
*                              Global Main Funcs                              *
*******************************************************************************/


void global_help(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    kwip SUBCOMMAND [options]\n");
    fprintf(stream, "\n");
    fprintf(stream, "Where subcommand is one of:\n\n");
    fprintf(stream, "  count:     Count kmers into sketches\n");
    fprintf(stream, "  kernel:    Calculate pairwise similarity between sketches\n");
    fprintf(stream, "\n");
    fprintf(stream, "To obtain help for a subcommand, use \"kwip SUBCOMMAND --help\"\n");
}

int main(int argc, char *argv[])
{
    if (argc < 2 || \
            strcmp(argv[1], "-h") == 0 || \
            strcmp(argv[1], "--help") == 0) {
        global_help(stdout);
        return EXIT_SUCCESS;
    }
    if (strcmp(argv[1], "-V") == 0 || \
            strcmp(argv[1], "--version") == 0) {
        kwip_print_version(stdout);
        return EXIT_SUCCESS;
    }

    if (strcmp(argv[1], "count") == 0) {
        return count_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "kernel") == 0) {
        return kernel_main(argc - 1, argv + 1);
    } else {
        global_help(stderr);
        return EXIT_FAILURE;
    }
}
