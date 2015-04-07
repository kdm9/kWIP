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

#include "d2thresh.hh"

namespace kmerclust
{
namespace metrics
{

DistanceCalcD2Thresh::
DistanceCalcD2Thresh()
{
    _threshold = 0;
}

float
DistanceCalcD2Thresh::
distance(khmer::CountingHash &a, khmer::CountingHash &b)
{
    size_t tab;
    size_t bin;
    std::vector<float> tab_scores;
    khmer::Byte **a_counts = a.get_raw_tables();
    khmer::Byte **b_counts = b.get_raw_tables();

    _check_hash_dimensions(a, b);

    //for (tab = 0; tab < a.n_tables(); tab++) {
    for (tab = 0; tab < 1; tab++) {
        uint64_t tab_dist = 0;
        khmer::Byte *A = a_counts[tab];
        khmer::Byte *B = b_counts[tab];
        for (bin = 0; bin < a.get_tablesizes()[tab]; bin++) {
            if (A[bin] > 0 && B[bin] > 0) {
                tab_dist += A[bin] * B[bin];
            }
        }
        tab_scores.push_back((float)tab_dist);
    }
    return tab_scores[0];
}

void
DistanceCalcD2Thresh::
set_threshold(unsigned int threshold)
{
    _threshold = threshold;
}

int
d2thresh_main(int argc, const char **argv)
{

    return EXIT_SUCCESS;
}


}} // end namespace kmerclust::metrics
