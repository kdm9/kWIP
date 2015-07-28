/*
 * ============================================================================
 *
 *       Filename:  entvecsums.cc
 *    Description:  Sum the inner product of a sample and the entropy vector.
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include <kernels/wip.hh>
#include <getopt.h>

class EntVecIPSummer: public kwip::metrics::WIPKernel
{

public:
    std::ostream               *tabstream = &std::cout;

    void
    calculate_pairwise          (std::vector<std::string> &hash_fnames);

    std::tuple<double, double>
    sample_entvec_sum           (khmer::CountingHash &sample);

};

void
EntVecIPSummer::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    num_samples = hash_fnames.size();
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < num_samples; i++) {
        add_hashtable(hash_fnames[i]);
        if (verbosity > 0) {
            #pragma omp critical
            {
                *outstream << "Loaded " << hash_fnames[i] << std::endl;
            }
        }
    }

    if (verbosity > 0) {
        *outstream << "Finished loading!" << std::endl;
        *outstream << "FPR: " << this->fpr() << std::endl;
    }

    _bin_entropies.clear();
    for (size_t tab = 0; tab < _n_tables; tab++) {
        _bin_entropies.emplace_back(_tablesizes[tab], 0.0);
        for (size_t bin = 0; bin < _tablesizes[tab]; bin++) {
            unsigned int bin_n_samples = _pop_counts[tab][bin];
            if (bin_n_samples == 0 || bin_n_samples == num_samples) {
                // Kmer not found in the population, or in all samples.
                // entropy will be 0, so bail out here
                _bin_entropies[tab][bin] = 0.0;
            } else {
                const float pop_freq = (float)bin_n_samples / (float)num_samples;
                _bin_entropies[tab][bin] = \
                        (pop_freq * -log2(pop_freq)) +
                        ((1 - pop_freq) * -log2(1 - pop_freq));
            }
        }
    }

    if (sample_names.empty()) {
        for (size_t i = 0; i < num_samples; i++) {
            size_t idx = hash_fnames[i].find_last_of("/");
            std::string base;
            if (idx != std::string::npos) {
                base = hash_fnames[i].substr(idx + 1);
            } else {
                base = hash_fnames[i];
            }
            base = base.substr(0, base.find_first_of("."));
            sample_names.push_back(std::string(base));
        }
    }

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < num_samples; i++) {
        kwip::CountingHashShrPtr ht = _get_hash(hash_fnames[i]);

        double sumentsamp = 0;
        double sumsamp = 0;
        std::tie(sumsamp, sumentsamp) = sample_entvec_sum(*ht);
        #pragma omp critical
        {
            *tabstream << sample_names[i] << "\t" << sumentsamp << "\t"
                       << sumsamp << "\t" << std::endl;
        }
    }
    if (verbosity > 0) {
        *outstream << "Done all!" << std::endl;
    }
}

std::tuple<double, double>
EntVecIPSummer::
sample_entvec_sum(khmer::CountingHash &sample)
{
    khmer::Byte **counts = sample.get_raw_tables();
    std::vector<double>     entvecsums;
    double count_sum = 0.0;

    for (size_t tab = 0; tab < _n_tables; tab++) {
        double countentvec_sum = 0.0;
        if (count_sum == 0.0) {
            // Count doesn't change between tables, hopefully
            for (size_t bin = 0; bin < _tablesizes[tab]; bin++) {
                count_sum += counts[tab][bin];
            }
        }
        for (size_t bin = 0; bin < _tablesizes[tab]; bin++) {
            float bin_entropy = _bin_entropies[tab][bin];
            float freq = counts[tab][bin] / count_sum;
            countentvec_sum += freq * bin_entropy;
        }
        entvecsums.push_back(countentvec_sum);
    }
    return std::make_tuple(count_sum, kwip::vec_min(entvecsums));
}


static std::string prog_name = "kwip-entvec";
static std::string cli_opts = "t:o:hVvq";

static const struct option cli_long_opts[] = {
    { "threads",    required_argument,  NULL,   't' },
    { "tabout",     required_argument,  NULL,   'o' },
    { "help",       no_argument,        NULL,   'h' },
    { "version",    no_argument,        NULL,   'V' },
    { "verbose",    no_argument,        NULL,   'v' },
    { "quiet",      no_argument,        NULL,   'q' },
};

static std::vector<std::string>
cli_help {
"-t, --threads       Number of threads to utilise. [default N_CPUS]",
"-o, --tabout        Output for tab-delimited sum table. [default stdout]",
"-h, --help          Print this help message.",
"-V, --version       Print the version string.",
"-v, --verbose       Increase verbosity. May or may not acutally do anything.",
"-q, --quiet         Execute silently but for errors.",
};

void
print_cli_help()
{
    using namespace std;

    cerr << endl;

    cerr << "USAGE: " << prog_name << " [options] hashes" << endl << endl;
    cerr << "OPTIONS:"<< endl;
    for (const auto &str: cli_help) {
        cerr << str << endl;
    }
    cerr << endl
         << "Each sample's oxli Countgraph should be specified after arguments:"
         << endl
         << prog_name << " [options] sample1.ct sample2.ct ... sampleN.ct"
         << endl;
}

int
main(int argc, char *argv[])
{
    EntVecIPSummer              evs;
    int                         option_idx      = 0;
    int                         c               = 0;
    std::string                 tabfile_name    = "";
    std::ofstream               tabfile;
    std::vector<std::string>    filenames;

    while ((c = getopt_long(argc, argv, cli_opts.c_str(), cli_long_opts,
                            &option_idx)) > 0) {
        switch (c) {
            case 't':
                evs.num_threads = atol(optarg);
                break;
            case 'o':
                tabfile_name = optarg;
                break;
            case 'v':
                evs.verbosity = 2;
                break;
            case 'q':
                evs.verbosity = 0;
                break;
            case 'h':
                print_cli_help();
                return EXIT_SUCCESS;
            case 'V':
                kwip::print_version();
                return EXIT_SUCCESS;
                break;
            case '?':
                print_cli_help();
                return EXIT_FAILURE;
        }
    }

    // Ensure we have at least two counting hashes to work with
    if (optind + 1 >= argc) {
        print_cli_help();
        return EXIT_FAILURE;
    }

    for (int i = optind; i < argc; i++) {
        filenames.push_back(std::string(argv[i]));
    }

    // Save matrices to files if we've been asked to
    if (tabfile_name.size() > 0) {
        tabfile.open(tabfile_name);
        evs.tabstream = &tabfile;
    }

    // Do the pairwise distance calculation
    evs.calculate_pairwise(filenames);

    tabfile.close();
    return EXIT_SUCCESS;
}
