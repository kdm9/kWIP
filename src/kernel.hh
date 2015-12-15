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

#ifndef KERNEL_HH
#define KERNEL_HH


#include <cmath>
#include <cassert>
#include <memory>
#include <limits>
#include <iostream>
#include <string>


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

#include <oxli/counting.hh> // liboxli countgraphs

#include "kwip-utils.hh"
#include "lrucache.hpp"


namespace kwip
{

typedef std::shared_ptr<khmer::CountingHash> CountingHashShrPtr;
typedef cache::lru_cache<std::string, CountingHashShrPtr> CountingHashCache;

class Kernel
{
protected:
    float                     **_kernel_mat;
    float                     **_distance_mat;
    MatrixXd                    _kernel_m;
    MatrixXd                    _normkernel_m;
    MatrixXd                    _distance_mt;
    int                         _num_threads;
    omp_lock_t                  _kernel_mat_lock;
    omp_lock_t                  _distance_mat_lock;
    CountingHashCache           _hash_cache;
    omp_lock_t                  _hash_cache_lock;

    // Ensure `a` and `b` have the same counting hash dimensions. Throws an
    // exception if they are not.
    virtual void
    _check_hash_dimensions     (khmer::CountingHash        &a,
                                khmer::CountingHash        &b);
    CountingHashShrPtr
    _get_hash                  (std::string                &filename);

    void
    _make_matrices             ();

    void
    _print_mat                 (std::ostream               &outstream,
                                float                     **matrix);


public:
    int                         verbosity;
    size_t                      num_samples;
    std::vector<std::string>    sample_names;
    const std::string           name = "Base Class";
    const std::string           blurb = "A generic base class for kernels.";
    std::ostream               *outstream = &std::cerr;

    Kernel                      ();
    ~Kernel                     ();

    // Calculate the kernel between two counting hashes
    virtual float
    kernel                      (khmer::CountingHash   &a,
                                 khmer::CountingHash   &b);

    // Calculate the kernel between all pairs of counting hashes in parallel
    virtual void
    calculate_pairwise          (std::vector<std::string> &hash_fnames);

    virtual void
    print_kernel_mat            (std::ostream          &outstream=std::cout);

    virtual void
    print_distance_mat          (std::ostream          &outstream=std::cout);

    virtual void
    kernel_to_distance          ();

    virtual float **
    get_kernel_matrix           ();

    virtual float **
    get_distance_matrix         ();

    virtual void
    load                        (std::istream          &instream);

    virtual void
    save                        (std::ostream          &outstream);

    void
    set_num_threads             (int                    num_threads);

};

} // end namespace kwip

#endif /* KERNEL_HH */
