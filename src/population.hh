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

#ifndef POPULATION_HH
#define POPULATION_HH


#include "kernel.hh"

namespace kmerclust
{

template<typename bin_tp>
class KernelPopulation : public Kernel
{

protected:
    bin_tp **_pop_counts;
    size_t _n_tables;
    std::vector<khmer::HashIntoType> _tablesizes;
    std::vector<uint64_t> _table_sums;
    omp_lock_t _pop_table_lock;

    void
    _check_pop_counts          (khmer::CountingHash        &ht);
public:
    KernelPopulation();

    ~KernelPopulation();

#if 0
    virtual void
    save                       (std::string                 filename);

    virtual void
    load                       (std::string                 filename);
#endif

    virtual void
    add_hashtable              (std::string                &hash_fname);

    virtual void
    calculate_pairwise         (std::vector<std::string>   &hash_fnames);

    double
    fpr                        ();
};

} // end namespace kmerclust

#endif /* POPULATION_HH */
