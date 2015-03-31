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


#include <klust.hh>

#include "lrucache-refcnt.hh"

#include <vector>
#include <queue>

#include <utility>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <omp.h>
#include <string.h>

using namespace std;
using namespace khmer;
using namespace kmerclust;
using namespace kmerclust::metrics;
using namespace refcounted_lru_cache;

static int hits = 0;
static int misses = 0;

std::shared_ptr<CountingHash>
get_hash(lru_cache<const char *, std::shared_ptr<CountingHash>> &cache, const char *path)
{
    std::shared_ptr<CountingHash> ret;
    while (1) {
        try {
            ret = cache.get(path);
            __sync_fetch_and_add(&hits, 1);
            return ret;
        } catch (std::range_error &err) {
            __sync_fetch_and_add(&misses, 1);
            std::shared_ptr<CountingHash> ht = std::make_shared<CountingHash>(1, 1);
            ht->load(path);
            cache.put(path, ht);

        }
    }
}

template<typename DistMeasure>
int
run_main (int argc, const char *argv[])
{
    DistMeasure distcalc;
    map<pair<int, int>, double> distances;


    CountingHash ht(1, 1);
    ht.load(argv[1]);
    distcalc.add_hashtable(ht);
    cerr << "Loaded " << argv[1] << endl;

    #pragma omp parallel for shared(distcalc)
    for (int i = 2; i < argc; i++) {
        CountingHash ht(1, 1);
    	ht.load(argv[i]);
        distcalc.add_hashtable(ht);
        #pragma omp critical
        {
            cerr << "Loaded " << argv[i] << endl;
        }
    }

    cerr << "Finished loading!" << endl;
    cerr << "FPR: " << distcalc.fpr() << endl;
}

int
main (int argc, char *argv[])
{
    if (argc < 4) {
        cerr << "USAGE: " << argv[0] << " " << argv[1] << \
            " <hashtable> ..." << endl;
        return EXIT_FAILURE;
    }

    DistanceCalcD2 distcalc;
    std::vector<std::string> filenames;

    for (int i = 1; i < argc; i++) {
        filenames.push_back(std::string(argv[i]));
    }

#if 0
    CountingHash ht(1, 1);
    ht.load(argv[1]);
    distcalc.add_hashtable(ht);
    cerr << "Loaded " << argv[1] << endl;
    #pragma omp parallel for shared(distcalc)
    for (int i = 2; i < argc; i++) {
        CountingHash ht(1, 1);
    	ht.load(argv[i]);
        distcalc.add_hashtable(ht);
        #pragma omp critical
        {
            cerr << "Loaded " << argv[i] << endl;
        }
    }

    cerr << "Finished loading!" << endl;
    cerr << "FPR: " << distcalc.fpr() << endl;
#endif
    distcalc.calculate_pairwise(filenames);
    distcalc.print_dist_mat();
    return EXIT_SUCCESS;
}
