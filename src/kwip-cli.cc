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

#include "kwip.hh"

#include <getopt.h>

using namespace kmerclust;  // Imports the base classes
using namespace kmerclust::metrics;  // Imports the KernelXXX classes

static std::string prog_name = "kwip";
static std::string cli_opts = "t:k:d:hUVvq";

static const struct option cli_long_opts[] = {
    { "threads",    required_argument,  NULL,   't' },
    { "kernel",     required_argument,  NULL,   'k' },
    { "distance",   required_argument,  NULL,   'd' },
    { "help",       no_argument,        NULL,   'h' },
    { "unweighted", no_argument,        NULL,   'U' },
    { "version",    no_argument,        NULL,   'V' },
    { "verbose",    no_argument,        NULL,   'v' },
    { "quiet",      no_argument,        NULL,   'q' },
};

static std::vector<std::string>
cli_help {
" -t, --threads     Number of threads to utilise [default N_CPUS].",
" -k, --kernel      Output file for the kernel matrix. [default None]",
" -d, --distance    Output file for the distance matrix. [default stdout]",
" -U, --unweighted  Use the unweighted inner proudct kernel. [default off]",
" -V, --version     Print the version string.",
" -v, --verbose     Increase verbosity. May or may not acutally do anything.",
" -q, --quiet       Execute silently but for errors.",
};

void
print_cli_help()
{
    using namespace std;

    print_version();
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

    while ((c = getopt_long(argc, argv, cli_opts.c_str(), cli_long_opts,
                            &option_idx)) > 0) {
        switch (c) {
            case 't':
                kernel.num_threads = atol(optarg);
                break;
            case 'k':
                kern_out_name = optarg;
                break;
            case 'd':
                dist_out_name = optarg;
                break;
            case 'v':
                kernel.verbosity = 2;
                break;
            case 'q':
                kernel.verbosity = 0;
                break;
            // This section is for the global options
            case 'h':
            case 'V':
            case 'U':
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
    if (dist_out_name.size() > 0) {
        dist_out.open(dist_out_name);
    }
    if (kern_out_name.size() > 0) {
        kern_out.open(kern_out_name);
    }

    // Do the pairwise distance calculation
    kernel.calculate_pairwise(filenames);

    // Only save the kernel distance if we have been given a file, or -
    if (kern_out_name == "-") {
        kernel.print_kernel_mat();
    } else if (kern_out_name.size() > 0) {
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

int
main (int argc, char *argv[])
{
    bool unweighted = false;
    prog_name = std::string(argv[0]);

    if (argc < 2) {
        print_cli_help();
        return EXIT_FAILURE;
    }

    int c;
    while ((c = getopt_long(argc, argv, cli_opts.c_str(), cli_long_opts,
                            NULL)) > 0) {
        switch (c) {
            case 'h':
                print_cli_help();
                return EXIT_SUCCESS;
            case 'V':
                print_version();
                return EXIT_SUCCESS;
            case 'U':
                unweighted = true;
                break;
            // This section is the kernel options, which we parse in the kernel
            // main function above.
            case 't':
            case 'k':
            case 'd':
            case 'q':
            case 'v':
                break;
            case '?':
                print_cli_help();
                return EXIT_FAILURE;
        }
    }

    if (unweighted) {
        return run_main<KernelD2>(argc, argv);
    }
    return run_main<KernelD2Ent>(argc, argv);
}
