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

namespace metrics
{
class KernelD2Thresh;
}

class Kernel
{
protected:
    int                         _n_threads;
    int                         _verbosity;
    size_t                      _n_samples;
    float                     **_kernel_mat;
    float                     **_distance_mat;
    omp_lock_t                  _kernel_mat_lock;
    omp_lock_t                  _distance_mat_lock;
    std::vector<std::string>    _sample_names;
    CountingHashCache           _hash_cache;
    omp_lock_t                  _hash_cache_lock;
    std::string                 _kernel_name = "Base class";

    // Ensure `a` and `b` have the same counting hash dimensions. Throws an
    // exception if they are not.
    virtual void
    _check_hash_dimensions     (khmer::CountingHash        &a,
                                khmer::CountingHash        &b);
    CountingHashShrPtr
    _get_hash                  (std::string                &filename);

    void
    _make_matrices             (size_t                      n_samples);

    void
    _print_mat                 (std::ostream               &outstream,
                                float                     **matrix);


public:
    Kernel                     ();
    ~Kernel                    ();

    // Calculate the kernel between two counting hashes
    virtual float
    kernel                     (khmer::CountingHash        &a,
                                khmer::CountingHash        &b);

    // Calculate the kernel between all pairs of counting hashes in parallel
    virtual void
    calculate_pairwise         (std::vector<std::string>   &hash_fnames);

    // Sets names of each sample in the hash
    virtual void
    set_sample_names           (std::vector<std::string>   &sample_names);

    virtual void
    set_num_threads            (int                         n_threads);

    virtual void
    print_kernel_mat           (std::ostream               &outstream);
    virtual void
    print_kernel_mat           ();

    virtual void
    print_distance_mat         (std::ostream               &outstream);
    virtual void
    print_distance_mat         ();

    virtual void
    kernel_to_distance         ();

    void
    set_n_threads              (int                         n_threads);

    void
    set_verbosity              (int                         verbosity);

    const std::string           blurb = "Base class";
};

} // end namespace kmerclust

#endif /* DISTANCE_HH */
