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

#include "d2ent.hh"

namespace kmerclust
{
namespace metrics
{

void
KernelD2Ent::
add_hashtable(const std::string &hash_fname)
{
    khmer::CountingHash ht(1, 1);
    khmer::Byte **counts;

    ht.load(hash_fname);
    _check_pop_counts(ht);

    counts = ht.get_raw_tables();

    for (size_t i = 0; i < 1; i++) {
        uint64_t tab_count = 0;
        // Save these here to avoid dereferencing twice below.
        uint16_t *this_popcount = _pop_counts[i];
        khmer::Byte *this_count = counts[i];
        for (size_t j = 0; j < _tablesizes[i]; j++) {
            if (this_count[j] > 0) {
                __sync_fetch_and_add(&(this_popcount[j]), 1);
            }
            tab_count += this_count[j];
        }
        __sync_fetch_and_add(&_table_sums[i], tab_count);
    }
}
void
KernelD2Ent::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    num_samples = hash_fnames.size();
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < num_samples; i++) {
        add_hashtable(hash_fnames[i]);
        if (verbosity > 0) {
            #pragma omp critical
            {
                std::cerr << "Loaded " << hash_fnames[i] << std::endl;
            }
        }
    }

    if (verbosity > 0) {
        std::cerr << "Finished loading!" << std::endl;
        std::cerr << "FPR: " << this->fpr() << std::endl;
    }

    double sum_bin_entropy = 0.0;
    _bin_entropies.clear();
    _bin_entropies.assign(_tablesizes[0], 0.0);
    for (size_t bin = 0; bin < _tablesizes[0]; bin++) {
        unsigned int bin_n_samples = _pop_counts[0][bin];
        if (bin_n_samples == 0 || bin_n_samples == num_samples) {
            // Kmer not found in the population, or in all samples.
            // entropy will be 0, so bail out here
            _bin_entropies[bin] = 0.0;
        } else {
            const float pop_freq = (float)bin_n_samples / (float)num_samples;
            _bin_entropies[bin] = (pop_freq * -log2(pop_freq)) +
                                  ((1 - pop_freq) * -log2(1 - pop_freq));
        }
        sum_bin_entropy += _bin_entropies[bin];
    }
    if (verbosity > 0) {
        std::cerr << "Finished loading!" << std::endl;
        std::cerr << "FPR: " << this->fpr() << std::endl;
    }
    // Do the kernel calculation per Kernel's implementation
    Kernel::calculate_pairwise(hash_fnames);
}

float
KernelD2Ent::
kernel(khmer::CountingHash &a, khmer::CountingHash &b)
{
    std::vector<float> tab_kernels;
    khmer::Byte **a_counts = a.get_raw_tables();
    khmer::Byte **b_counts = b.get_raw_tables();

    _check_hash_dimensions(a, b);

    for (size_t tab = 0; tab < 1; tab++) {
        float tab_kernel = 0.0;
        double sum_a = 0, sum_b = 0;
        for (size_t bin = 0; bin < _tablesizes[tab]; bin++) {
            sum_a += a_counts[tab][bin];
            sum_b += b_counts[tab][bin];
        }
        for (size_t bin = 0; bin < _tablesizes[tab]; bin++) {
            float bin_entropy = _bin_entropies[bin];
            float a_freq = a_counts[tab][bin] / sum_a;
            float b_freq = b_counts[tab][bin] / sum_b;
            tab_kernel += a_freq * b_freq * bin_entropy;
        }
        tab_kernels.push_back(tab_kernel);
    }
    return tab_kernels[0];
}

}} // end namespace kmerclust::metrics
