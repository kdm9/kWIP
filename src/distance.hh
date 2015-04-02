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

#ifndef DISTANCE_HH
#define DISTANCE_HH

#include <cmath>
#include <memory>

#ifdef _OPENMP
    #include <omp.h>
#else
    #warning "All the data and no openmp makes Jack a dull boy"
    #define omp_lock_t int
    #define omp_init_lock (void)
    #define omp_set_lock (void)
    #define omp_unset_lock (void)
    #define omp_destroy_lock (void)
    #define omp_get_max_threads(x) (1)
#endif

#include <counting.hh>

#include "lrucache.hpp"

namespace kmerclust
{

typedef std::shared_ptr<khmer::CountingHash> CountingHashShrPtr;
typedef cache::lru_cache<std::string, CountingHashShrPtr> CountingHashCache;

class DistanceCalc
{

protected:
    int _n_threads;
    size_t _n_samples;
    float **_dist_mat;
    omp_lock_t _dist_mat_lock;
    std::vector<std::string> _sample_names;
    CountingHashCache _hash_cache;
    omp_lock_t _hash_cache_lock;

    // Ensure `a` and `b` have the same counting hash dimensions. Throws an
    // exception if they are not.
    virtual void
    _check_hash_dimensions     (khmer::CountingHash        &a,
                                khmer::CountingHash        &b);
    CountingHashShrPtr
    _get_hash                  (std::string                &filename);

public:
    DistanceCalc               ();
    ~DistanceCalc              ();

    // Caclulate the distance between two counting hashes
    virtual float
    distance                   (khmer::CountingHash        &a,
                                khmer::CountingHash        &b);

    virtual void
    calculate_pairwise         (std::vector<std::string>   &hash_fnames);

    virtual void
    set_sample_names           (std::vector<std::string>   &sample_names);

    virtual void
    set_num_threads            (int                         n_threads);

    virtual void
    print_dist_mat             (std::ostream               &outstream);

    virtual void
    print_dist_mat             ()
    { print_dist_mat(std::cout); }

};

} // end namespace kmerclust

#endif /* DISTANCE_HH */
