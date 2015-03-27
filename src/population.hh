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


#include "distance.hh"

#ifndef POPULATION_HH
#define POPULATION_HH

namespace kmerclust
{

class DistanceCalcPopulation : public DistanceCalc
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

    void _make_tables(std::vector<khmer::HashIntoType> &tablesizes);

public:
    DistanceCalcPopulation()
    {
        _pop_counts = NULL;
        _n_tables = 0;
    }

    ~DistanceCalcPopulation();

    virtual void save(std::string filename) { }

    virtual void load(std::string filename) { }

    void add_hashtable(khmer::CountingHash &ht);

    double fpr();
};

} // end namespace kmerclust

#endif /* POPULATION_HH */
