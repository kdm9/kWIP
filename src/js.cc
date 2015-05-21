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

#include "js.hh"

namespace kmerclust
{
namespace metrics
{

float
KernelJS::
kernel(khmer::CountingHash &a, khmer::CountingHash &b)
{
    std::vector<float> tab_scores;
    khmer::Byte **a_counts = a.get_raw_tables();
    khmer::Byte **b_counts = b.get_raw_tables();
    std::vector<khmer::HashIntoType> tablesizes = a.get_tablesizes();

    _check_hash_dimensions(a, b);

    //for (size_t tab = 0; tab < a.n_tables(); tab++) {
    for (size_t tab = 0; tab < 1; tab++) {
        float small_log_offset = 0.0;
        uint64_t a_sum = 0;
        uint64_t b_sum = 0;
        float tab_kernel = 0;
        khmer::Byte *A = a_counts[tab];
        khmer::Byte *B = b_counts[tab];

        // Calculate the sums absolute values of each table
        for (size_t bin = 0; bin < tablesizes[tab]; bin++) {
            // NB: the bins of a CountingHash can never be negative, so there's
            // no need to use `abs()`
            a_sum += A[bin];
            b_sum += B[bin];
        }
        const uint64_t max_sum = a_sum > b_sum ? a_sum : b_sum;
        small_log_offset = 1.0 / (float)max_sum;

        for (size_t bin = 0; bin < tablesizes[tab]; bin++) {
            const float a_freq = (float)A[bin] / (float)a_sum;
            const float b_freq = (float)B[bin] / (float)b_sum;

            const float total = a_freq + b_freq + small_log_offset;

            const float a_log = a_freq * log2((a_freq + small_log_offset)/total);
            const float b_log = b_freq * log2((b_freq + small_log_offset)/total);

            tab_kernel += -(a_log + b_log);
        }
        tab_scores.push_back((float)tab_kernel);
    }
    return tab_scores[0];
}


}} // end namespace kmerclust::metrics
