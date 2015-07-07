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

#ifndef D2ENT_HH
#define D2ENT_HH


#include "population.hh"

#include <sstream>
#include <iostream>
#include <limits>

namespace kwip
{
namespace metrics
{

class KernelD2Ent : public KernelPopulation<uint16_t>
{

public:
    void
    add_hashtable               (const std::string     &hash_fname);

    float
    kernel                      (khmer::CountingHash   &a,
                                 khmer::CountingHash   &b);

    void
    calculate_pairwise          (std::vector<std::string> &hash_fnames);

    void
    calculate_entropy_vector    (std::vector<std::string> &hash_fnames);

    const std::string       blurb =
            "D2Ent Kernel\n"
            "\n"
            "This kernel calculates the D2 Entropy kernel, the inner product\n"
            "bin frequencies (bin value / total number of kmers) weighted by\n"
            "the Shannon entropy of presence/absence over all samples (i.e.\n"
            "the proportion of all samples with a non-zero count of a bin).\n"
            "\n"
            "D2ent = sum (A_f * B_f * P_ent)\n"
            "Where:\n"
            "A_f, B_f = A[i] / sum A\n"
            "P_ent = (P[i] / sum P) * - log2(P[i] / sum P)\n"
            " where P is a vector of bin occurance across all samples\n";

    void
    load                        (std::istream       &instream);

    void
    save                        (std::ostream       &outstream);

private:
    std::vector<float>      _bin_entropies;
    const std::string       _file_sig="kWIP_BinEntVector";
};

}} // end namespace kwip::metrics


#endif /* D2ENT_HH */
