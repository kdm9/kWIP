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

using namespace kwip;  // Imports the base classes
using namespace kwip::metrics;  // Imports the KernelXXX classes

static std::string prog_name = "kwip";
static std::string cli_opts = "t:k:d:w:hCUVvq";

static const struct option cli_long_opts[] = {
    { "threads",    required_argument,  NULL,   't' },
    { "kernel",     required_argument,  NULL,   'k' },
    { "distance",   required_argument,  NULL,   'd' },
    { "weights",    required_argument,  NULL,   'w' },
    { "help",       no_argument,        NULL,   'h' },
    { "calc-weights", no_argument,      NULL,   'C' },
    { "unweighted", no_argument,        NULL,   'U' },
    { "version",    no_argument,        NULL,   'V' },
    { "verbose",    no_argument,        NULL,   'v' },
    { "quiet",      no_argument,        NULL,   'q' },
};

static std::vector<std::string>
cli_help {
"-t, --threads       Number of threads to utilise. [default N_CPUS]",
"-k, --kernel        Output file for the kernel matrix. [default None]",
"-d, --distance      Output file for the distance matrix. [default stdout]",
"-U, --unweighted    Use the unweighted inner proudct kernel. [default off]",
"-w, --weights       Bin weight vector file (input, or output w/ -C).",
"-C, --calc-weights  Calculate only the bin weight vector, not kernel matrix.",
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

template<typename KernelImpl>
int
run_pwcalc(int argc, char *argv[])
{
    KernelImpl                  kernel;
    int                         option_idx      = 0;
    int                         c               = 0;
    std::string                 dist_out_name   = "";
    std::string                 kern_out_name   = "";
    std::ofstream               dist_out;
    std::ofstream               kern_out;
    std::ifstream               weights_file;
    std::vector<std::string>    filenames;

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
            case 'w':
                weights_file.open(optarg);
                break;
            // This section is for the global options
            case 'h':
            case 'V':
            case 'U':
            case 'C':
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

    if (weights_file.is_open()) {
        kernel.load(weights_file);
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
run_precalc(int argc, char *argv[])
{
    KernelD2Ent                 kernel;
    int                         option_idx      = 0;
    int                         c               = 0;
    std::ofstream               weights_file;
    std::string                 weights_file_name;
    std::vector<std::string>    filenames;

    while ((c = getopt_long(argc, argv, cli_opts.c_str(), cli_long_opts,
                            &option_idx)) > 0) {
        switch (c) {
            case 't':
                kernel.num_threads = atol(optarg);
                break;
            case 'v':
                kernel.verbosity = 2;
                break;
            case 'q':
                kernel.verbosity = 0;
                break;
            case 'w':
                weights_file_name = optarg;
                weights_file.open(optarg);
                if (!weights_file.is_open()) {
                    std::cerr << "Error opening weights file '"
                              << weights_file_name << "'" << std::endl;
                    print_cli_help();
                    return EXIT_FAILURE;
                }
                break;
            // This section is for the pairwise calculation main options
            case 'k':
            case 'd':
            // This section is for the global options
            case 'h':
            case 'V':
            case 'U':
            case 'C':
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

    if (weights_file_name.size() < 1) {
        std::cerr << "Weight vector file must be supplied with '-C'"
                  << std::endl;
        print_cli_help();
        return EXIT_FAILURE;
    }

    for (int i = optind; i < argc; i++) {
        filenames.push_back(std::string(argv[i]));
    }
    kernel.calculate_entropy_vector(filenames);

    kernel.save(weights_file);
    return EXIT_SUCCESS;
}

int
main (int argc, char *argv[])
{
    bool unweighted = false;
    bool precalc_weights = false;
    int retval = 0;
    int opt = 0;
    prog_name = std::string(argv[0]);

    if (argc < 2) {
        print_cli_help();
        return EXIT_FAILURE;
    }

    int c;
    while ((c = getopt_long(argc, argv, cli_opts.c_str(), cli_long_opts,
                            &opt)) > 0) {
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
            case 'C':
                precalc_weights = true;
                break;
            // This section is the kernel options, which we parse in the kernel
            // main function above.
            case 't':
            case 'k':
            case 'd':
            case 'q':
            case 'v':
            case 'w':
                break;
            case '?':
                print_version();
                print_cli_help();
                return EXIT_FAILURE;
        }
    }
    // Reset so that getops works later on.
    optind = 0;

    if (precalc_weights) {
        retval = run_precalc(argc, argv);
    } else if (unweighted) {
        retval = run_pwcalc<KernelD2>(argc, argv);
    } else {
        retval = run_pwcalc<KernelD2Ent>(argc, argv);
    }
    return retval;
}
