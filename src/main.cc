//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#include "khmer.hh"
#include "counting.hh"
#include "distance.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashbits.hh"
#include "labelhash.hh"
#include "khmer_exception.hh"

#include "lrucache.hh"

#include <vector>
#include <queue>

#include <thread>
#include <mutex>
#include <condition_variable>

#include <omp.h>
#include <string.h>

using namespace std;
using namespace khmer;
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

int
main (int argc, const char *argv[])
{
    CountingHashDistanceCalcD2pop distcalc;
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

#if 0
    for (int i = 1; i < argc; i++) {
        CountingHash ht1(1, 1);
        CountingHash ht2(1, 1);
        ht1.load(argv[i]);
        for (int j = 1; j < argc; j++) {
            double dist = 0.0;
            const pair<ssize_t, ssize_t> ij(i, j);
            if (i > j) {
                continue;
            }
            ht2.load(argv[j]);
            num_loaded++;
            dist = distcalc.distance(ht1, ht2);
            distances[ij] = dist;
	    cerr << i << " x " << j << " done!" << endl;
        }
    }
#else
    lru_cache<const char *, std::shared_ptr<CountingHash>> cache(omp_get_num_procs() * 4);

    #pragma omp parallel for shared(distcalc, distances, cache) schedule(dynamic)
    for (int i = 1; i < argc; i++) {
        std::shared_ptr<CountingHash> ht1;
        std::shared_ptr<CountingHash> ht2;
        #pragma omp critical
        {
            ht1 = get_hash(cache, argv[i]);
        }
        for (int j = 1; j < argc; j++) {
            double dist = 0.0;
            pair<int, int> ij(i, j);
            if (i > j) {
                // Skip the bottom half of the distance matrix.
                continue;
            }
            #pragma omp critical
            {
                ht2 = get_hash(cache, argv[j]);
            }
            dist = distcalc.distance(*ht1, *ht2);
            #pragma omp critical
            {
                distances[ij] = dist;
                cerr << i << " x " << j << endl;
                cache.unget(argv[j]);
            }
        }
        #pragma omp critical
        {
            cache.unget(argv[i]);
        }
    }
#endif

    for (ssize_t i = 1; i < argc; i++) {
        cout << "\t" << basename(argv[i]);
    }
    cout << endl;
    for (ssize_t i = 1; i < argc; i++) {
        cout << basename(argv[i]);
        for (ssize_t j = 1; j < argc; j++) {
            if (i > j) {
                pair<ssize_t, ssize_t> ij(j, i);
                cout << "\t" << distances[ij];
            } else {
                pair<ssize_t, ssize_t> ij(i, j);
                cout << "\t" << distances[ij];
            }
        }
        cout << endl;
    }
    cerr << "Done! " << hits << " hits, " << misses << " misses in cache." << endl;
    return EXIT_SUCCESS;
} /* ----------  end of function main  ---------- */
