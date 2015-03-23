//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#ifndef DISTANCE_HH
#define DISTANCE_HH

#include "counting.hh"
#include <cmath>

namespace khmer
{

class CountingHashDistanceCalc
{

protected:
    void _check_hash_dimensions(CountingHash &a, CountingHash &b)
    {
        size_t i;
        bool ok = true;
	std::vector<HashIntoType> a_tsz = a.get_tablesizes(),
				  b_tsz = b.get_tablesizes();

        if (a.ksize() != b.ksize()) ok = false;
        if (a.n_tables() != b.n_tables()) ok = false;
        for (i = 0; i <a.n_tables(); i++) {
            if (a_tsz[i] != b_tsz[i]) ok = false;
        }
        if (!ok) {
            throw khmer_exception("Hash dimensions and k-size not equal");
        }
    }

public:
    virtual float distance(CountingHash &a, CountingHash &b)
    {
	    return 0.0;
    }

};

class CountingHashDistanceCalcD2 : public CountingHashDistanceCalc
{
public:
    float distance(CountingHash &a, CountingHash &b)
    {
        size_t tab;
        size_t bin;
        std::vector<float> tab_scores;
	khmer::Byte **a_counts = a.get_raw_tables(),
	            **b_counts = b.get_raw_tables();

        _check_hash_dimensions(a, b);

        for (tab = 0; tab < a.n_tables(); tab++) {
            float tab_dist = 0.0;
            for (bin = 0; bin < a.get_tablesizes()[tab]; bin++) {
                tab_dist += a_counts[tab][bin] * b_counts[tab][bin];
            }
            tab_scores.push_back(tab_dist);
        }
        return tab_scores[0];
    }
};

class CountingHashDistanceCalcPopulation : public CountingHashDistanceCalc
{

protected:
    uint16_t **_pop_counts;
    size_t _n_tables;
    std::vector<HashIntoType> _tablesizes;
    std::vector<uint64_t> _table_sums;

    bool _have_tables()
    {
        return (_pop_counts != NULL);
    }

    void _make_tables(std::vector<HashIntoType> &tablesizes)
    {
        size_t i;

        _tablesizes = tablesizes;
        _n_tables = tablesizes.size();
        _pop_counts = new uint16_t*[_n_tables];
        for (i = 0; i < _n_tables; i++) {
            _pop_counts[i] = new uint16_t[tablesizes[i]];
            memset(_pop_counts[i], 0, tablesizes[i] * sizeof(uint16_t));
	    _table_sums.push_back(0);
        }
    }

public:
    CountingHashDistanceCalcPopulation()
    {
        _pop_counts = NULL;
	_n_tables = 0;
    }

    virtual ~CountingHashDistanceCalcPopulation()
    {
        size_t i;

	if (_pop_counts == NULL) return;
        for (i = 0; i < _n_tables; i++) {
            delete[] _pop_counts[i];
        }
        delete[] _pop_counts;
    }

    virtual void save(std::string filename)
    {
    }

    virtual void load(std::string filename)
    {
    }

    void add_hashtable(CountingHash &ht)
    {
        size_t i, j;
        khmer::Byte **counts = ht.get_raw_tables();

        if (!_have_tables()) {
            std::vector<HashIntoType> tablesizes = ht.get_tablesizes();
            _make_tables(tablesizes);
        }
        for (i = 0; i < _n_tables; i++) {
            uint64_t tab_count = 0;
            for (j = 0; j < _tablesizes[i]; j++) {
                __sync_fetch_and_add(&_pop_counts[i][j], counts[i][j]);
                tab_count += counts[i][j];
            }
            __sync_fetch_and_add(&_table_sums[i], tab_count);
        }
    }
};

class CountingHashDistanceCalcD2pop : public CountingHashDistanceCalcPopulation
{

public:

    float distance(CountingHash &a, CountingHash &b)
    {
        size_t tab;
        size_t bin;
        std::vector<float> tab_dists;
	khmer::Byte **a_counts = a.get_raw_tables(),
	            **b_counts = b.get_raw_tables();

        _check_hash_dimensions(a, b);

        for (tab = 0; tab < 1; tab++) {
            float tab_dist = 0.0;
            uint64_t sum_a = 0, sum_b = 0;

            for (bin = 0; bin < _tablesizes[tab]; bin++) {
                sum_a += a_counts[tab][bin];
                sum_b += b_counts[tab][bin];
            }
            for (bin = 0; bin < _tablesizes[tab]; bin++) {
                float cent_a, cent_b; // Centered word counts
                float bin_freq; // Population freq of bin

                bin_freq =  _pop_counts[tab][bin] / (float)_table_sums[tab];

                cent_a = a_counts[tab][bin] - (sum_a * bin_freq);
                cent_b = b_counts[tab][bin] - (sum_b * bin_freq);

                tab_dist += cent_a * cent_b;
            }
            tab_dists.push_back(tab_dist);
        }

        return tab_dists[0];
    }
};

}

#endif /* DISTANCE_HH */
