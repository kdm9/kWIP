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

namespace kmerclust
{
namespace metrics
{

class KernelD2Ent : public KernelPopulation<uint16_t>
{

public:
    KernelD2Ent                 ();
    ~KernelD2Ent                ();

    void
    add_hashtable               (const std::string     &hash_fname);

    float kernel                (khmer::CountingHash   &a,
                                 khmer::CountingHash   &b);
private:
    std::vector<float>      _bin_entropies;
    omp_lock_t              _bin_entropy_vec_lock;
};

}} // end namespace kmerclust::metrics


#endif /* D2ENT_HH */
