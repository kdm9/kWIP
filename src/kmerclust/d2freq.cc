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

#include "d2freq.hh"

namespace kmerclust
{
namespace metrics
{

float
KernelD2freq::kernel(khmer::CountingHash &a, khmer::CountingHash &b)
{
    std::vector<float> tab_kernels;
    khmer::Byte **a_counts = a.get_raw_tables();
    khmer::Byte **b_counts = b.get_raw_tables();
    std::vector<khmer::HashIntoType> tablesizes = a.get_tablesizes();

    _check_hash_dimensions(a, b);

    for (size_t tab = 0; tab < 1; tab++) {
        uint64_t tabsz = tablesizes[tab];
        float tab_kernel = 0;
        double sum_a = 0, sum_b = 0;
        khmer::Byte *A = a_counts[tab];
        khmer::Byte *B = b_counts[tab];
        for (size_t bin = 0; bin < tabsz; bin++) {
            sum_a += A[bin];
            sum_b += B[bin];
        }
        for (size_t bin = 0; bin < tabsz; bin++) {
            float a_freq = A[bin] / sum_a;
            float b_freq = B[bin] / sum_b;
            tab_kernel += a_freq * b_freq;
        }
        tab_kernels.push_back(tab_kernel);
    }
    return tab_kernels[0];
}

}} // end namespace kmerclust::metrics
