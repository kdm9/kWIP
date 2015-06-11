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
    _kernel_mat(NULL),
    _distance_mat(NULL),
    _hash_cache(1),
    verbosity(1),
    num_samples(0)
{
    omp_init_lock(&_kernel_mat_lock);
    omp_init_lock(&_distance_mat_lock);
    omp_init_lock(&_hash_cache_lock);
    num_threads = omp_get_max_threads();
    _hash_cache = CountingHashCache(num_threads*2 + 1);
}

Kernel::
~Kernel()
{
    // Lock these to avoid a race. Very unlikely but doesn't hurt
    omp_set_lock(&_kernel_mat_lock);
    omp_set_lock(&_distance_mat_lock);
    omp_set_lock(&_hash_cache_lock);
    // Free the kernel matrix if it exists
    if (_kernel_mat != NULL) {
        for (size_t i = 0; i < num_samples; i++) {
            delete[] _kernel_mat[i];
        }
        delete[] _kernel_mat;
    }
    // Free dist matrix if it exists
    if (_distance_mat != NULL) {
        for (size_t i = 0; i < num_samples; i++) {
            delete[] _distance_mat[i];
        }
        delete[] _distance_mat;
    }
    // Destroy all locks
    omp_destroy_lock(&_kernel_mat_lock);
    omp_destroy_lock(&_distance_mat_lock);
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
_make_matrices()
{
    // Create the kernel matrix
    omp_set_lock(&_kernel_mat_lock);
    omp_set_lock(&_distance_mat_lock);
    if (_kernel_mat != NULL) {
        for (size_t i = 0; i < num_samples; i++) {
            delete[] _kernel_mat[i];
        }
        delete[] _kernel_mat;
    }
    _kernel_mat = new float *[num_samples];
    for (size_t i = 0; i < num_samples; i++) {
        _kernel_mat[i] = new float[num_samples];
    }

    // Create the distance matrix
    if (_distance_mat != NULL) {
        for (size_t i = 0; i < num_samples; i++) {
            delete[] _distance_mat[i];
        }
        delete[] _distance_mat;
    }
    _distance_mat = new float *[num_samples];
    for (size_t i = 0; i < num_samples; i++) {
        _distance_mat[i] = new float[num_samples];
    }
    omp_unset_lock(&_kernel_mat_lock);
    omp_unset_lock(&_distance_mat_lock);
}

void
Kernel::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    num_samples = hash_fnames.size();

    _make_matrices();

    if (sample_names.empty()) {
        for (size_t i = 0; i < num_samples; i++) {
            char *fname = strdup(hash_fnames[i].c_str());
            assert(fname != NULL);
            std::string base{basename(fname)};
            base = base.substr(0, base.find_first_of("."));
            sample_names.push_back(std::string(base));
            free(fname);
        }
    }

    #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
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
            _kernel_mat[i][j] = kernel;
            _kernel_mat[j][i] = kernel;
            if (verbosity > 0) {
                #pragma omp critical
                {
                    std::cerr << i << " x " << j << " done!" << std::endl;
                }
            }
        }
    }
    if (verbosity > 0) {
        std::cerr << "Done all!" << std::endl;
    }
}

void
Kernel::
_print_mat(std::ostream &outstream, float **matrix)
{
    // Avoid a segfault
    if (matrix == NULL) {
        throw std::runtime_error("Invalid matrix provided");
    }
    for (size_t i = 0; i < num_samples; i++) {
        if (matrix[i] == NULL) {
            throw std::runtime_error("Invalid matrix provided");
        }
    }

    // Use numerals as indices if there's no names provided
    if (sample_names.empty()) {
        for (size_t i = 0; i < num_samples; i++) {
            // This cast is required to fix an error with Intel compilers
            sample_names.push_back(std::to_string((unsigned long long)i));
        }
    }

    // Header row
    outstream << "."; // Marker for the far top left cell
    for (const auto &sample: sample_names) {
        outstream << "\t" << sample;
    }
    outstream << std::endl;

    // The matrix itself
    for (size_t i = 0; i < num_samples; i++) {
        outstream << sample_names[i];
        for (size_t j = 0; j < num_samples; j++) {
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
print_distance_mat(std::ostream &outstream)
{
    kernel_to_distance();
    _print_mat(outstream, _distance_mat);
}

void
Kernel::
kernel_to_distance()
{
    std::vector<float> diag(num_samples);

    if (_distance_mat == NULL) {
        throw std::runtime_error("No distance matrix exists");
    }
    if (_kernel_mat == NULL) {
        throw std::runtime_error("No kernel matrix exists");
    }

    float **tmp_mat = new float *[num_samples];
    for (size_t i = 0; i < num_samples; i++) {
        tmp_mat[i] = new float[num_samples];
    }

    // Store the diagonal of the matrix
    for (size_t i = 0; i < num_samples; i++) {
        diag[i] = _kernel_mat[i][i];
    }

    for (size_t i = 0; i < num_samples; i++) {
        for (size_t j = 0; j < num_samples; j++) {
            float this_val = _kernel_mat[i][j];
            float norm_factor = sqrt(diag[i] * diag[j]);
            tmp_mat[i][j] = this_val / norm_factor;
        }
    }

    for (size_t i = 0; i < num_samples; i++) {
        for (size_t j = 0; j < num_samples; j++) {
            float dist = tmp_mat[i][i] + tmp_mat[j][j] - 2 * tmp_mat[i][j];
            if (dist > 0.0) {
                dist = sqrt(dist);
            }
            _distance_mat[i][j] = dist;
        }
    }

    // Free the temporary matrix
    for (size_t i = 0; i < num_samples; i++) {
        delete[] tmp_mat[i];
    }
    delete[] tmp_mat;
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
