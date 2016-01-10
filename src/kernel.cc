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

#include "kernel.hh"

#include <Eigen/Eigenvalues>

namespace kwip
{

Kernel::
Kernel() :
    _kernel_m(1,1),
    _hash_cache(1),
    verbosity(1),
    num_samples(0)
{
    omp_init_lock(&_hash_cache_lock);
    _num_threads = omp_get_max_threads();
    _hash_cache = CountingHashCache(_num_threads + 1);
}

Kernel::
~Kernel()
{
    // Lock these to avoid a race. Very unlikely but doesn't hurt
    omp_set_lock(&_hash_cache_lock);
    omp_destroy_lock(&_hash_cache_lock);
}

float
Kernel::
kernel(khmer::CountingHash &a, khmer::CountingHash &b)
{
    _check_hash_dimensions(a, b);
    return 0.0;
}


void
Kernel::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    num_samples = hash_fnames.size();

    _kernel_m.resize(num_samples, num_samples);
    _kernel_m.fill(0);

    if (sample_names.empty()) {
        for (size_t i = 0; i < num_samples; i++) {
            size_t idx = hash_fnames[i].find_last_of("/");
            std::string base;
            if (idx != std::string::npos) {
                base = hash_fnames[i].substr(idx + 1);
            } else {
                base = hash_fnames[i];
            }
            std::vector<std::string> exts {
                ".kh",
                ".ct",
                ".cg",
                ".countgraph",
            };
            size_t ext_idx = std::string::npos;
            for (const auto &ext: exts) {
                ext_idx = base.find(ext);
                if (ext_idx != std::string::npos) {
                    base.erase(ext_idx);
                    break;
                }
            }
            sample_names.push_back(std::string(base));
        }
    }

    #pragma omp parallel for schedule(dynamic) num_threads(_num_threads)
    for (size_t i = 0; i < num_samples; i++) {
        CountingHashShrPtr ht1 = _get_hash(hash_fnames[i]);
        for (size_t j = 0; j < num_samples; j++) {
            float kernel = 0.0;
            if (i > j) {
                // Skip calculating the bottom half of the matrix
                continue;
            }
            CountingHashShrPtr ht2 = _get_hash(hash_fnames[j]);
            kernel = this->kernel(*ht1, *ht2);
            // Fill in both halves of the matrix
            _kernel_m(i, j) = kernel;
            _kernel_m(j, i) = kernel;
            if (verbosity > 0) {
                #pragma omp critical
                {
                    *outstream << i + 1 << " x " << j + 1 << " done!"
                               << std::endl;
                }
            }
        }
    }
    if (verbosity > 0) {
        *outstream << "Done all!" << std::endl;
    }
}

void
Kernel::
print_kernel_mat(std::ostream &outstream)
{
    std::cerr << "pkm ksz " << _kernel_m.size() << "\n";
    if (_kernel_m.size() == 1) {
        // No kernel has been calculated
        throw std::runtime_error("No kernel matrix exists");
    }
    print_lsmat(_kernel_m, outstream, sample_names);
}

void
Kernel::
get_kernel_matrix(MatrixXd &mat)
{
    if (_kernel_m.size() == 1) {
        // No kernel has been calculated
        throw std::runtime_error("No kernel matrix exists");
    }
    mat = _kernel_m;
}

void
Kernel::
get_norm_kernel_matrix(MatrixXd &mat)
{
    if (_kernel_m.size() == 1) {
        // No kernel has been calculated
        throw std::runtime_error("No kernel matrix exists");
    }
    normalise_matrix(mat, _kernel_m);
}

void
Kernel::
get_distance_matrix(MatrixXd &mat)
{
    if (_kernel_m.size() == 1) {
        // No kernel has been calculated
        throw std::runtime_error("No kernel matrix exists");
    }
    MatrixXd dist(num_samples, num_samples);
    kernel_to_distance(mat, _kernel_m);
}

void
Kernel::
print_distance_mat(std::ostream &outstream)
{
    if (_kernel_m.size() == 1) {
        // No kernel has been calculated
        throw std::runtime_error("No kernel matrix exists");
    }
    MatrixXd dist;
    get_distance_matrix(dist);
    print_lsmat(dist, outstream, sample_names);
}


void
Kernel::
load(std::istream &instream)
{
    (void)instream;
}

void
Kernel::
save(std::ostream &outstream)
{
    (void)outstream;
}

void
Kernel::
set_num_threads(int num_threads)
{
    _num_threads = num_threads;

    // Update the hash cache size
    _hash_cache = CountingHashCache(_num_threads + 1);
}

CountingHashShrPtr
Kernel::
_get_hash(std::string &filename)
{
    CountingHashShrPtr ret;
    omp_set_lock(&_hash_cache_lock);
    while (1) {
        try {
            ret = _hash_cache.get(filename);
            omp_unset_lock(&_hash_cache_lock);
            return ret;
        } catch (std::range_error &err) {
            CountingHashShrPtr ht = \
                    std::make_shared<khmer::CountingHash>(1, 1);
            khmer::CountingHashFile::load(filename, *ht);
            _hash_cache.put(filename, ht);
        }
    }
}

void
Kernel::
_check_hash_dimensions(khmer::CountingHash &a, khmer::CountingHash &b)
{
    size_t i;
    bool ok = true;
    std::vector<khmer::HashIntoType> a_tsz = a.get_tablesizes();
    std::vector<khmer::HashIntoType> b_tsz = b.get_tablesizes();

    if (a.ksize() != b.ksize()) {
        ok = false;
    }
    if (a.n_tables() != b.n_tables()) {
        ok = false;
    }
    for (i = 0; i <a.n_tables(); i++) {
        if (a_tsz[i] != b_tsz[i])  {
            ok = false;
        }
    }
    if (!ok) {
        throw std::runtime_error("Hash dimensions and k-size not equal");
    }
}

} // namespace kwip
