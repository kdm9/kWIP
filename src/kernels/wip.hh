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

class WIPKernel : public KernelPopulation<uint16_t>
{

public:
    void
    add_hashtable               (const std::string     &hash_fname);

    float
    kernel                      (const khmer::CountingHash   &a,
                                 const khmer::CountingHash   &b);

    void
    calculate_pairwise          (std::vector<std::string> &hash_fnames);

    void
    calculate_entropy_vector    (std::vector<std::string> &hash_fnames);

    const std::string       blurb =
            "WIP Kernel\n"
            "\n"
            "This kernel calculates the Shannon Entropy Weighted Inner Product\n"
            "kernel, the inner product of bin frequencies (bin value / total\n"
            "number of kmers) weighted by the Shannon entropy of presence-\n"
            "absence over all samples (i.e. the proportion of all samples with\n"
            "a non-zero count of a bin).\n";

    void
    load                        (std::istream       &instream);

    void
    save                        (std::ostream       &outstream);

protected:
    std::vector<std::vector<float>>  _bin_entropies;
    const std::string       _file_sig="kWIP_BinEntVector";
};

}} // end namespace kwip::metrics


#endif /* D2ENT_HH */
