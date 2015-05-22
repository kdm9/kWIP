/*
 * Copyright 2015 Kevin Murray <spam@kdmurray.id.au>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kmerclust.hh"

#include <getopt.h>

using namespace kmerclust;  // Imports the base classes
using namespace kmerclust::metrics;  // Imports the KernelXXX classes

static std::string cli_opts = "t:k:d:hVvq";

static std::vector<struct option>
cli_long_opts = {
    { "threads",    required_argument,  NULL,   't' },
    { "kernel",     required_argument,  NULL,   'k' },
    { "distance",   required_argument,  NULL,   'd' },
    { "help",       no_argument,        NULL,   'h' },
    { "version",    no_argument,        NULL,   'V' },
    { "verbose",    no_argument,        NULL,   'v' },
    { "quiet",      no_argument,        NULL,   'q' },
};

static std::vector<std::string>
cli_help = {
    "  -t, --threads    Number of threads to utilise [default N_CPUS].",
    "  -k, --kernel     Output file for the kernel matrix. [default None]",
    "  -d, --distance   Output file for the distance matrix. [default stdout]",
    "  -V, --version    Print the version string.",
    "  -v, --verbose    Increase verbosity.",
    "  -q, --quiet      Execute silently but for errors.",
};

void
print_version()
{
    using namespace std;

    cerr << "kmerclust version " << KMERCLUST_VERSION << endl;
}

void
print_cli_help(std::string prog, std::string kernel)
{
    using namespace std;

    print_version();
    cerr << endl;

    cerr << "USAGE: " << prog << " " << kernel << " [options] hashes" << endl
         << endl;
    cerr << "OPTIONS:"<< endl;
    for (auto str: cli_help) {
        cerr << str << endl;
    }
    cerr << endl;
    cerr << "Khmer CountingHashes should be specified as files after arguments, e.g:" << endl;
    cerr << prog << " " << kernel << " *.kh" << endl;
}

template<typename KernelImpl>
int
run_main(int argc, char *argv[])
{
    KernelImpl                  kernel;
    int                         option_idx      = 0;
    int                         c               = 0;
    std::string                 dist_out_name   = "";
    std::string                 kern_out_name   = "";
    std::ofstream               dist_out;
    std::ofstream               kern_out;
    std::vector<std::string>    filenames;
    std::string                 prog            = argv[0];
    std::string                 kernel_abbrev   = argv[1];
    const char                 *opts;
    struct option              *long_opts;
    size_t                      n_options       = cli_long_opts.size();

    // Kernel-specific CLI options
    if (std::is_same<KernelImpl, KernelD2Thresh>()) {
        cli_opts += "T:";
        cli_long_opts.push_back({ "threshold", required_argument, NULL, 'T' });
        cli_help.push_back(
            "  -T, --threshold  Threshold for the inner product calculation.");
    }

    // remove argv[0], which is just the binary name. This makes getopt think
    // the kernel name is the binary name, and makes everything work nicely.
    argc--;
    argv++;

    // Create option C string
    opts = cli_opts.c_str();

    // Create long options array
    long_opts = new struct option[n_options + 1];
    for (size_t i = 0; i < n_options; i++) {
        long_opts[i] = cli_long_opts[i];
    }
    long_opts[n_options] = {NULL, 0, NULL, 0}; // Add sentinel

    while ((c = getopt_long(argc, argv, opts, long_opts, &option_idx)) > 0) {
        switch (c) {
            case 't':
                kernel.set_n_threads(atol(optarg));
                break;
            case 'k':
                kern_out_name = optarg;
                break;
            case 'd':
                dist_out_name = optarg;
                break;
            case 'h':
                print_cli_help(prog, kernel_abbrev);
                std::cerr << std::endl << "---- Kernel Details ----"
                          << std::endl;
                std::cerr << std::endl << kernel.blurb << std::endl;
                delete [] long_opts;
                return EXIT_SUCCESS;
            case 'V':
                print_version();
                delete [] long_opts;
                return EXIT_SUCCESS;
            case 'v':
                kernel.set_verbosity(2);
                break;
            case 'q':
                kernel.set_verbosity(0);
                break;
            case 'T':
                if (std::is_same<KernelImpl, KernelD2Thresh>()) {
                    // Cast to a d2pop
                    // Because fuck C++, that's why
                    Kernel         *base  = \
                                    dynamic_cast<Kernel *>(&kernel);
                    KernelD2Thresh *d2pop = \
                                    dynamic_cast<KernelD2Thresh *>(base);
                    d2pop->set_threshold(atoi(optarg));
                } else {
                    // It's an error if any other kernel has '-T'
                    print_cli_help(prog, kernel_abbrev);
                    delete [] long_opts;
                    return EXIT_FAILURE;
                }
            case '?':
                // Getopt long prints its own error msg
                print_cli_help(prog, kernel_abbrev);
                delete [] long_opts;
                return EXIT_FAILURE;
        }
    }
    delete [] long_opts;

    // Ensure we have at least two counting hashes to work with
    if (optind + 1 >= argc) {
        print_cli_help(prog, kernel_abbrev);
    }

    for (int i = optind; i < argc; i++) {
        filenames.push_back(std::string(argv[i]));
    }

    // Save matrices to files if we've been asked to
    if (dist_out_name.size() > 0) {
        dist_out.open(dist_out_name);
    }
    if (kern_out_name.size() > 0) {
        kern_out.open(kern_out_name);
    }

    // Do the pairwise distance calculation
    kernel.calculate_pairwise(filenames);

    // Only save the kernel distance if we have been given a file, or -
    if (dist_out_name == "-") {
        kernel.print_kernel_mat();
    } else if (dist_out_name.size() > 0) {
        kernel.print_kernel_mat(kern_out);
    }
    // Always save the distance matrix, to stdout if we don't have a file
    if (dist_out_name.size() > 0 && dist_out_name != "-") {
        kernel.print_distance_mat(dist_out);
    } else {
        kernel.print_distance_mat();
    }

    kern_out.close();
    dist_out.close();
    return EXIT_SUCCESS;
}

void
print_valid_kernels()
{
    std::cerr << "Valid kernels are:" << std::endl;
    std::cerr << "  d2" << std::endl
              << "  d2pop" << std::endl
              << "  d2ent" << std::endl
              << "  d2thresh" << std::endl
              << "  d2freq" << std::endl
              << "  js" << std::endl;
}

int
main (int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "USAGE: " << argv[0]
                  << " <kernel> [options] <hashtable> ..." << std::endl;
        print_valid_kernels();
        std::cerr << "See " << argv[0] << " <kernel> -h for further help."
                  << std::endl;
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "d2") == 0) {
        return run_main<KernelD2>(argc, argv);
    } else if (strcmp(argv[1], "d2pop") == 0) {
        return run_main<KernelD2pop>(argc, argv);
    } else if (strcmp(argv[1], "d2ent") == 0) {
        return run_main<KernelD2Ent>(argc, argv);
    } else if (strcmp(argv[1], "d2thresh") == 0) {
        return run_main<KernelD2Thresh>(argc, argv);
    } else if (strcmp(argv[1], "d2freq") == 0) {
        return run_main<KernelD2freq>(argc, argv);
    } else if (strcmp(argv[1], "js") == 0) {
        return run_main<KernelJS>(argc, argv);
    }

    // If we get to here, we have an error
    std::cerr << "ERROR: Invalid kernel name " << argv[1] << std::endl
              << std::endl;
    print_valid_kernels();
    return EXIT_FAILURE;
}
