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

DistanceCalcPopulation::
DistanceCalcPopulation():
    _pop_counts(NULL),
    _n_tables(0)
{
    omp_init_lock(&_pop_table_lock);
}

DistanceCalcPopulation::
~DistanceCalcPopulation()
{
    omp_destroy_lock(&_pop_table_lock);
    if (_pop_counts != NULL) {
        for (size_t i = 0; i < _n_tables; i++) {
            delete[] _pop_counts[i];
        }
        delete[] _pop_counts;
    }
}

void
DistanceCalcPopulation::
add_hashtable(khmer::CountingHash &ht)
{
    khmer::Byte **counts = ht.get_raw_tables();

    omp_set_lock(&_pop_table_lock);
    if (!_have_tables()) {
        std::vector<khmer::HashIntoType> tablesizes = ht.get_tablesizes();
        _tablesizes = tablesizes;
        _n_tables = tablesizes.size();
        _pop_counts = new uint16_t*[_n_tables];
        for (size_t i = 0; i < _n_tables; i++) {
            _pop_counts[i] = new uint16_t[tablesizes[i]];
            memset(_pop_counts[i], 0, tablesizes[i] * sizeof(uint16_t));
        _table_sums.push_back(0);
        }
    }
    omp_unset_lock(&_pop_table_lock);
    // Below here is threadsafe, I think

    for (size_t i = 0; i < _n_tables; i++) {
        uint64_t tab_count = 0;
        // Save these here to avoid dereferencing twice below.
        uint16_t *this_popcount = _pop_counts[i];
        khmer::Byte *this_count = counts[i];
        for (size_t j = 0; j < _tablesizes[i]; j++) {
            __sync_fetch_and_add(&this_popcount[j], this_count[j]);
            tab_count += this_count[j];
        }
        __sync_fetch_and_add(&_table_sums[i], tab_count);
    }
}

bool
DistanceCalcPopulation::
_have_tables()
{
    return (_pop_counts != NULL);
}

void
DistanceCalcPopulation::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    _n_samples = hash_fnames.size();

    #pragma omp parallel for num_threads(_n_threads)
    for (size_t i = 0; i < _n_samples; i++) {
        khmer::CountingHash ht(1, 1);
        ht.load(hash_fnames[i]);
        add_hashtable(ht);
        #pragma omp critical
        {
            std::cerr << "Loaded " << hash_fnames[i] << std::endl;
        }
    }

    std::cerr << "Finished loading!" << std::endl;
    std::cerr << "FPR: " << this->fpr() << std::endl;

    // Do the distance calculation per DistanceCalc's implementation
    DistanceCalc::calculate_pairwise(hash_fnames);
}

double
DistanceCalcPopulation::
fpr()
{
    double fpr = 1;
    std::vector<double> tab_counts(_n_tables, 0);

    for (size_t i = 0; i < _n_tables; i++) {
        uint64_t tab_count = 0;
        #pragma omp parallel for num_threads(_n_threads)
        for (size_t j = 0; j < _tablesizes[i]; j++) {
            tab_count += _pop_counts[i][j] > 0 ? 1 : 0;
        }
        tab_counts[i] = (double)tab_count / _tablesizes[i];
    }

    for (size_t i = 0; i < _n_tables; i++) {
        fpr *= tab_counts[i];
    }
    return fpr;
}

} // end namespace kmerclust


