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

#include "population.hh"

namespace kmerclust
{

DistanceCalcPopulation::
~DistanceCalcPopulation()
{
    size_t i;

    if (_pop_counts == NULL) return;
    for (i = 0; i < _n_tables; i++) {
        delete[] _pop_counts[i];
    }
    delete[] _pop_counts;
}

void
DistanceCalcPopulation::
_make_tables(std::vector<khmer::HashIntoType> &tablesizes)
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

void
DistanceCalcPopulation::
add_hashtable(khmer::CountingHash &ht)
{
    size_t i, j;
    khmer::Byte **counts = ht.get_raw_tables();

    if (!_have_tables()) {
        std::vector<khmer::HashIntoType> tablesizes = ht.get_tablesizes();
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

double
DistanceCalcPopulation::
fpr()
{
    size_t i, j;
    double fpr = 1;
    std::vector<double> tab_counts(_n_tables, 0);

    #pragma omp parallel for num_threads(_n_threads)
    for (i = 0; i < _n_tables; i++) {
        uint64_t tab_count = 0;
        for (j = 0; j < _tablesizes[i]; j++) {
            tab_count += _pop_counts[i][j] > 0 ? 1 : 0;
        }
        tab_counts[i] = (double)tab_count / _tablesizes[i];
    }

    for (i = 0; i < _n_tables; i++) {
        fpr *= tab_counts[i];
    }
    return fpr;
}

} // end namespace kmerclust


