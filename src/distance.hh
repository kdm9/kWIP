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

#ifndef DISTANCE_HH
#define DISTANCE_HH

#include "counting.hh"
#include <cmath>

namespace kmerclust
{

class CountingHashDistanceCalc
{

protected:
    void _check_hash_dimensions(khmer::CountingHash &a, khmer::CountingHash &b);

public:
    virtual float distance(khmer::CountingHash &a, khmer::CountingHash &b)
    {
	return 0.0;
    }

};

class CountingHashDistanceCalcD2 : public CountingHashDistanceCalc
{
public:
    float distance(khmer::CountingHash &a, khmer::CountingHash &b);
};

class CountingHashDistanceCalcPopulation : public CountingHashDistanceCalc
{

protected:
    uint16_t **_pop_counts;
    size_t _n_tables;
    std::vector<khmer::HashIntoType> _tablesizes;
    std::vector<uint64_t> _table_sums;

    bool _have_tables()
    {
        return (_pop_counts != NULL);
    }

    void _make_tables(std::vector<khmer::HashIntoType> &tablesizes)
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

    virtual void save(std::string filename) { }

    virtual void load(std::string filename) { }

    void add_hashtable(khmer::CountingHash &ht);
};

class CountingHashDistanceCalcD2pop : public CountingHashDistanceCalcPopulation
{

public:

    float distance(khmer::CountingHash &a, khmer::CountingHash &b);
};

}

#endif /* DISTANCE_HH */
