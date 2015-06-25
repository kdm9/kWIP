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

template<typename bin_tp>
KernelPopulation<bin_tp>::
KernelPopulation():
    _pop_counts(NULL),
    _n_tables(0)
{
    omp_init_lock(&_pop_table_lock);
}

template<typename bin_tp>
KernelPopulation<bin_tp>::
~KernelPopulation()
{
    omp_destroy_lock(&_pop_table_lock);
    if (_pop_counts != NULL) {
        for (size_t i = 0; i < _n_tables; i++) {
            delete[] _pop_counts[i];
        }
        delete[] _pop_counts;
    }
}

template<typename bin_tp>
void
KernelPopulation<bin_tp>::
add_hashtable(const std::string &hash_fname)
{
    khmer::CountingHash ht(1, 1);
    khmer::Byte **counts;

    ht.load(hash_fname);
    _check_pop_counts(ht);

    counts = ht.get_raw_tables();

    for (size_t i = 0; i < _n_tables; i++) {
        uint64_t tab_count = 0;
        // Save these here to avoid dereferencing twice below.
        bin_tp *this_popcount = _pop_counts[i];
        khmer::Byte *this_count = counts[i];
        for (size_t j = 0; j < _tablesizes[i]; j++) {
            __sync_fetch_and_add(&(this_popcount[j]), this_count[j]);
            tab_count += this_count[j];
        }
        __sync_fetch_and_add(&_table_sums[i], tab_count);
    }
}

template<typename bin_tp>
void
KernelPopulation<bin_tp>::
_check_pop_counts(khmer::CountingHash &ht)
{
    omp_set_lock(&_pop_table_lock);
    if (_pop_counts == NULL) {
        std::vector<khmer::HashIntoType> tablesizes = ht.get_tablesizes();
        _tablesizes = tablesizes;
        _n_tables = tablesizes.size();
        _pop_counts = new bin_tp*[_n_tables];
        for (size_t i = 0; i < _n_tables; i++) {
            _pop_counts[i] = new bin_tp[tablesizes[i]];
            memset(_pop_counts[i], 0, tablesizes[i] * sizeof(bin_tp));
        _table_sums.push_back(0);
        }
    }
    omp_unset_lock(&_pop_table_lock);
}


template<typename bin_tp>
void
KernelPopulation<bin_tp>::
_free_pop_counts()
{
    omp_set_lock(&_pop_table_lock);
    if (_pop_counts != NULL) {
        for (size_t i = 0; i < _n_tables; i++) {
            delete[] _pop_counts[i];
        }
        delete[] _pop_counts;
    }
    omp_unset_lock(&_pop_table_lock);
}


template<typename bin_tp>
void
KernelPopulation<bin_tp>::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < hash_fnames.size(); i++) {
        add_hashtable(hash_fnames[i]);
        if (verbosity > 0) {
            #pragma omp critical
            {
                *outstream << "Loaded " << hash_fnames[i] << std::endl;
            }
        }
    }

    if (verbosity > 0) {
        *outstream << "Finished loading!" << std::endl;
        *outstream << "FPR: " << this->fpr() << std::endl;
    }

    // Do the kernel calculation per Kernel's implementation
    Kernel::calculate_pairwise(hash_fnames);
}

template<typename bin_tp>
double
KernelPopulation<bin_tp>::
fpr()
{
    double fpr = 1;
    std::vector<double> tab_counts(_n_tables, 0);

    for (size_t i = 0; i < _n_tables; i++) {
        uint64_t tab_count = 0;
        for (size_t j = 0; j < _tablesizes[i]; j++) {
            tab_count += _pop_counts[i][j] > 0 ? 1 : 0;
        }
        tab_counts[i] = (double)tab_count / (double)_tablesizes[i];
    }

    for (size_t i = 0; i < _n_tables; i++) {
        fpr *= tab_counts[i];
    }
    return fpr;
}

#if 0
template<typename bin_tp>
void
KernelPopulation<bin_tp>::
save(const std::string &filename)
{
}

template<typename bin_tp>
void
KernelPopulation<bin_tp>::
load(const std::string &filename)
{
}
#endif

// Explicit compilation of standard types
template class KernelPopulation<uint8_t>;
template class KernelPopulation<uint16_t>;
template class KernelPopulation<uint32_t>;
template class KernelPopulation<uint64_t>;

} // end namespace kmerclust


