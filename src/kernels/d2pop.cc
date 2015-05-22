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

#include "d2pop.hh"

namespace kmerclust
{
namespace metrics
{

float
KernelD2pop::kernel(khmer::CountingHash &a, khmer::CountingHash &b)
{
    size_t tab;
    size_t bin;
    std::vector<float> tab_kernels;
    khmer::Byte **a_counts = a.get_raw_tables();
    khmer::Byte **b_counts = b.get_raw_tables();

    _check_hash_dimensions(a, b);

    for (tab = 0; tab < 1; tab++) {
        float tab_kernel = 0.0;
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

            tab_kernel += cent_a * cent_b;
        }
        tab_kernels.push_back(tab_kernel);
    }

    return tab_kernels[0];
}

}} // end namespace kmerclust::metrics
