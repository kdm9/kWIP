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

namespace kmerclust
{

void
CountingHashDistanceCalc::_check_hash_dimensions(khmer::CountingHash &a, khmer::CountingHash &b)
{
    size_t i;
    bool ok = true;
    std::vector<khmer::HashIntoType> a_tsz = a.get_tablesizes();
    std::vector<khmer::HashIntoType> b_tsz = b.get_tablesizes();

    if (a.ksize() != b.ksize()) ok = false;
    if (a.n_tables() != b.n_tables()) ok = false;
    for (i = 0; i <a.n_tables(); i++) {
        if (a_tsz[i] != b_tsz[i]) ok = false;
    }
    if (!ok) {
        throw "Hash dimensions and k-size not equal";
    }
}

float
CountingHashDistanceCalcD2::distance(khmer::CountingHash &a, khmer::CountingHash &b)
{
    size_t tab;
    size_t bin;
    std::vector<float> tab_scores;
    khmer::Byte **a_counts = a.get_raw_tables();
    khmer::Byte **b_counts = b.get_raw_tables();

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

void
CountingHashDistanceCalcPopulation::add_hashtable(khmer::CountingHash &ht)
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
CountingHashDistanceCalcPopulation::fpr()
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

float
CountingHashDistanceCalcD2pop::distance(khmer::CountingHash &a, khmer::CountingHash &b)
{
    size_t tab;
    size_t bin;
    std::vector<float> tab_dists;
    khmer::Byte **a_counts = a.get_raw_tables();
    khmer::Byte **b_counts = b.get_raw_tables();

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

} // namespace kmerclust
