#include <oxli/counting.hh>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <getopt.h>

void
usage(FILE *stream)
{
    fprintf(stream, "oxlicap -- cap counts in oxli countgraphs\n");
    fprintf(stream, "\n");
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    oxlicap -c CAP COUNTFILE\n");
    fprintf(stream, "\n");
    fprintf(stream, "CAP must be between 1 and 255 inclusive.\n");
}

int
cap_file(const char *filename, int cap)
{
    using namespace khmer;
    using std::vector;
    using std::min;
    Byte **counts;
    size_t tables;
    vector<HashIntoType> tablesizes;

    CountingHash ht(1, 1);
    CountingHashFile::load(filename, ht);

    counts = ht.get_raw_tables();
    tables = ht.n_tables();

    std::cerr << "Capping " << filename << " to " << cap << "\n";
    for (size_t i = 0; i < tables; i++) {
        size_t tablesize = ht.get_tablesizes()[i];
        for (size_t j = 0; j < tablesize; j++) {
            counts[i][j] = min(cap, (int)counts[i][j]);
        }
    }

    std::string outfile (filename);
    outfile += std::string(".capped.gz");
    std::cerr << "Done capping. Saving to " << outfile << "\n";
    CountingHashFile::save(outfile.c_str(), ht);
    std::cerr << "All Done\n";
    return 0;
}

int
main(int argc, char *argv[])
{
    int cap = 0;

    int c;
    while ((c = getopt(argc, argv, "c:")) > 0) {
        switch (c) {
            case 'c':
                cap = atoi(optarg);
                if (cap < 1 || cap > 255) {
                    std::cerr << "ERROR: cap must be between 1 and 255 inclusive.\n";
                    return EXIT_FAILURE;
                }
                break;
            case '?':
                usage(stderr);
                return EXIT_FAILURE;
        }
    }

    if (optind > argc - 1) {
        usage(stdout);
        return EXIT_SUCCESS;
    }

    return cap_file(argv[optind], cap);
}
