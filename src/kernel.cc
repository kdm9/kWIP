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

namespace kmerclust
{

Kernel::
Kernel() :
    _verbosity(1),
    _n_samples(0),
    _kernel_mat(NULL),
    _distance_mat(NULL),
    _hash_cache(1)
{
    omp_init_lock(&_kernel_mat_lock);
    omp_init_lock(&_distance_mat_lock);
    omp_init_lock(&_hash_cache_lock);
    _n_threads = omp_get_max_threads();
    _hash_cache = CountingHashCache(_n_threads*2 + 1);
}

Kernel::
~Kernel()
{
    omp_set_lock(&_kernel_mat_lock);
    if (_kernel_mat == NULL) {
        omp_unset_lock(&_kernel_mat_lock);
        return;
    }
    for (size_t i = 0; i < _n_samples; i++) {
        delete[] _kernel_mat[i];
        delete[] _distance_mat[i];
    }
    delete[] _kernel_mat;
    omp_unset_lock(&_kernel_mat_lock);
    omp_destroy_lock(&_kernel_mat_lock);
}

float
Kernel::
kernel(khmer::CountingHash &a, khmer::CountingHash &b)
{
    (void) a;
    (void) b;
    return 0.0;
}

void
Kernel::
_make_matrices()
{
    // Create the kernel matrix
    omp_set_lock(&_kernel_mat_lock);
    if (_kernel_mat != NULL) {
        for (size_t i = 0; i < _n_samples; i++) {
            delete[] _kernel_mat[i];
        }
        delete[] _kernel_mat;
    }
    _kernel_mat = new float *[_n_samples];
    for (size_t i = 0; i < _n_samples; i++) {
        _kernel_mat[i] = new float[_n_samples];
    }
    omp_unset_lock(&_kernel_mat_lock);

    // Create the distance matrix
    omp_set_lock(&_distance_mat_lock);
    if (_distance_mat != NULL) {
        for (size_t i = 0; i < _n_samples; i++) {
            delete[] _distance_mat[i];
        }
        delete[] _distance_mat;
    }
    _distance_mat = new float *[_n_samples];
    for (size_t i = 0; i < _n_samples; i++) {
        _distance_mat[i] = new float[_n_samples];
    }
    omp_unset_lock(&_distance_mat_lock);
}

void
Kernel::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    _n_samples = hash_fnames.size();

    _make_matrices();

    if (_sample_names.empty()) {
        for (size_t i = 0; i < _n_samples; i++) {
            char *fname = strdup(hash_fnames[i].c_str());
            std::string base(basename(fname));
            base = base.substr(0, base.find_first_of("."));
            _sample_names.push_back(std::string(base));
            free(fname);
        }
    }

    #pragma omp parallel for schedule(dynamic) num_threads(_n_threads)
    for (size_t i = 0; i < _n_samples; i++) {
        CountingHashShrPtr ht1 = _get_hash(hash_fnames[i]);
        for (size_t j = 0; j < _n_samples; j++) {
            float kernel = 0.0;
            if (i > j) {
                // Skip calculating the bottom half of the matrix
                continue;
            }
            CountingHashShrPtr ht2 = _get_hash(hash_fnames[j]);
            kernel = this->kernel(*ht1, *ht2);
            // Fill in both halves of the matrix
            _kernel_mat[i][j] = kernel;
            _kernel_mat[j][i] = kernel;
            if (_verbosity > 0) {
                #pragma omp critical
                {
                    std::cerr << i << " x " << j << " done!" << std::endl;
                }
            }
        }
    }
    if (_verbosity > 0) {
        std::cerr << "Done all!" << std::endl;
    }
}

void
Kernel::
set_sample_names(std::vector<std::string> &sample_names)
{
    for (auto name: sample_names) {
        _sample_names.push_back(name);
    }
}

void
Kernel::
set_num_threads(int n_threads)
{
    _n_threads = n_threads;
}

void
Kernel::
_print_mat(std::ostream &outstream, float **matrix)
{
    // Avoid a segfault
    if (matrix == NULL) {
        throw std::runtime_error("Invalid matrix provided");
    }
    for (size_t i = 0; i < _n_samples; i++) {
        if (matrix[i] == NULL) {
            throw std::runtime_error("Invalid matrix provided");
        }
    }

    // Use numerals as indices if there's no names provided
    if (_sample_names.empty()) {
        for (size_t i = 0; i < _n_samples; i++) {
            // This cast is required to fix an error with Intel compilers
            _sample_names.push_back(std::to_string((unsigned long long)i));
        }
    }

    // Header row
    outstream << "."; // Marker for the far top left cell
    for (size_t i = 0; i < _n_samples; i++) {
        outstream << "\t" << _sample_names[i];
    }
    outstream << std::endl;
    // The matrix itself
    for (size_t i = 0; i < _n_samples; i++) {
        outstream << _sample_names[i];
        for (size_t j = 0; j < _n_samples; j++) {
            outstream << "\t" << matrix[i][j];
        }
        outstream << std::endl;
    }
}

void
Kernel::
print_kernel_mat(std::ostream &outstream)
{
    _print_mat(outstream, _kernel_mat);
}

void
Kernel::
print_kernel_mat()
{
    _print_mat(std::cout, _kernel_mat);
}

void
Kernel::
print_distance_mat(std::ostream &outstream)
{
    kernel_to_distance();
    _print_mat(outstream, _distance_mat);
}

void
Kernel::
print_distance_mat()
{
    kernel_to_distance();
    _print_mat(std::cout, _distance_mat);
}

void
Kernel::
kernel_to_distance()
{
    std::vector<float> diag(_n_samples);

    if (_distance_mat == NULL) {
        throw std::runtime_error("No distance matrix exists");
    }
    if (_kernel_mat == NULL) {
        throw std::runtime_error("No kernel matrix exists");
    }

    float **tmp_mat = new float *[_n_samples];
    for (size_t i = 0; i < _n_samples; i++) {
        tmp_mat[i] = new float[_n_samples];
    }

    // Store the diagonal of the matrix
    for (size_t i = 0; i < _n_samples; i++) {
        diag[i] = _kernel_mat[i][i];
    }

    for (size_t i = 0; i < _n_samples; i++) {
        for (size_t j = 0; j < _n_samples; j++) {
            float this_val = _kernel_mat[i][j];
            float norm_factor = sqrt(diag[i] * diag[j]);
            tmp_mat[i][j] = this_val / norm_factor;
        }
    }

    for (size_t i = 0; i < _n_samples; i++) {
        for (size_t j = 0; j < _n_samples; j++) {
            float dist = tmp_mat[i][i] + tmp_mat[j][j] - 2 * tmp_mat[i][j];
            if (dist > 0.0) {
                dist = sqrt(dist);
            }
            _distance_mat[i][j] = dist;
        }
    }

    // Free the temporay matrix
    for (size_t i = 0; i < _n_samples; i++) {
        delete [] tmp_mat[i];
    }
    delete [] tmp_mat;
}

void
Kernel::
set_n_threads(int n_threads)
{
    _n_threads = n_threads;
}

void
Kernel::
set_verbosity(int verbosity)
{
    _verbosity = verbosity;
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
            ht->load(filename);
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

} // namespace kmerclust
