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

#include "wip.hh"

namespace kwip
{
namespace metrics
{

void
WIPKernel::
add_hashtable(const std::string &hash_fname)
{
    khmer::CountingHash ht(1, 1);
    khmer::Byte **counts;

    khmer::CountingHashFile::load(hash_fname, ht);
    _check_pop_counts(ht);

    counts = ht.get_raw_tables();

    for (size_t tab = 0; tab < _n_tables; tab++) {
        uint64_t tab_count = 0;
        // Save these here to avoid dereferencing twice below.
        uint16_t *this_popcount = _pop_counts[tab];
        khmer::Byte *this_count = counts[tab];
        for (size_t j = 0; j < _tablesizes[tab]; j++) {
            if (this_count[j] > 0) {
                __sync_fetch_and_add(&(this_popcount[j]), 1);
            }
            tab_count += this_count[j];
        }
        __sync_fetch_and_add(&_table_sums[tab], tab_count);
    }
}


void
WIPKernel::
calculate_entropy_vector(std::vector<std::string> &hash_fnames)
{
    num_samples = hash_fnames.size();
    if (verbosity > 0) {
        *outstream << "Calculating entropy weighting vector:" << std::endl;
        *outstream << "  - Loading hashes into a population frequency vector:"
                   << std::endl;
    }

    #pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
    for (size_t i = 0; i < num_samples; i++) {
        add_hashtable(hash_fnames[i]);
        if (verbosity > 0) {
            #pragma omp critical
            {
                *outstream << "    - Loaded '" << hash_fnames[i]  << "' ("
                           << i + 1 << ")" << std::endl;
            }
        }
    }

    if (verbosity > 0) {
        *outstream << " - Finished loading hashes!" << std::endl;
        *outstream << " - Occupancy rate of population hash: "
                   << this->fpr() << std::endl;
    }

    _bin_entropies.clear();
    for (size_t tab = 0; tab < _n_tables; tab++) {
        _bin_entropies.emplace_back(_tablesizes[tab], 0.0);
        for (size_t bin = 0; bin < _tablesizes[tab]; bin++) {
            // Number of samples in popn with non-zero for this bin
            unsigned int pop_count = _pop_counts[tab][bin];
            if (0 < pop_count && pop_count < num_samples) {
                const float pop_freq = (float)pop_count / (float)num_samples;
                // Shannon entropy is
                // sum for all states p_state * -log_2(p_state)
                // We have two states, present & absent
                _bin_entropies[tab][bin] =  \
                        (pop_freq * -log2(pop_freq)) +
                        ((1 - pop_freq) * -log2(1 - pop_freq));
            }
        }
    }
}


void
WIPKernel::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    // Only load samples and calculate the bin entropy vector if we don't have
    // it already
    if (_bin_entropies.size() == 0) {
        calculate_entropy_vector(hash_fnames);
    } else {
        num_samples = hash_fnames.size();
    }

    // Do the kernel calculation per Kernel's implementation
    Kernel::calculate_pairwise(hash_fnames);
}

float
WIPKernel::
kernel(khmer::CountingHash &a, khmer::CountingHash &b)
{
    std::vector<float>          tab_kernels;
    khmer::Byte               **a_counts = a.get_raw_tables();
    khmer::Byte               **b_counts = b.get_raw_tables();

    _check_hash_dimensions(a, b);

    for (size_t tab = 0; tab < _n_tables; tab++) {
        double tab_kernel = 0.0;
        double norm_a = 0, norm_b = 0;
        for (size_t bin = 0; bin < _tablesizes[tab]; bin++) {
            norm_a += (uint32_t)a_counts[tab][bin] * (uint32_t)a_counts[tab][bin];
            norm_b += (uint32_t)b_counts[tab][bin] * (uint32_t)b_counts[tab][bin];
        }
        norm_a = sqrt(norm_a);
        norm_b = sqrt(norm_b);
        for (size_t bin = 0; bin < _tablesizes[tab]; bin++) {
            uint32_t a = a_counts[tab][bin];
            uint32_t b = b_counts[tab][bin];

            if (a == 0 || b == 0 ) {
                continue;
            }
            double bin_entropy = _bin_entropies[tab][bin];
            double a_freq = a / norm_a;
            double b_freq = b / norm_b;
            tab_kernel += a_freq * b_freq * bin_entropy;
        }
        tab_kernels.push_back(tab_kernel);
    }
    return vec_min(tab_kernels);
}

void
WIPKernel::
load(std::istream &instream)
{
    std::string filesig;
    int64_t     hashsize;

    instream >> filesig;
    instream >> hashsize;

    if (filesig != _file_sig) {
        std::runtime_error("Input is not a kWIP WIP bin entropy vector");
    }
    if (hashsize <= 0) {
        std::ostringstream msg;
        msg << "Invalid number of bins: " <<  hashsize;
        std::runtime_error(msg.str());
    }

    _bin_entropies.clear();
    _bin_entropies.emplace_back(hashsize, 0.0);
    for (ssize_t i = 0; i < hashsize; i++) {
        size_t idx;
        instream >> idx;
        instream >> _bin_entropies[0][i];
    }
}

void
WIPKernel::
save(std::ostream &outstream)
{
    if (_bin_entropies.size() < 1) {
        std::runtime_error("There is no bin entropy vector to save");
    }

    outstream.precision(std::numeric_limits<float>::digits10);
    outstream << _file_sig << "\t" << _bin_entropies[0].size() << "\n";
    for (size_t i = 0; i < _bin_entropies[0].size(); i++) {
        outstream << i << "\t" << _bin_entropies[0][i] << "\n";
    }
}

}} // end namespace kwip::metrics
