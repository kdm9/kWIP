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

#include "distance.hh"

namespace kmerclust
{

DistanceCalc::
DistanceCalc() :
    _n_samples(0),
    _dist_mat(NULL),
    _hash_cache(1)
{
    omp_init_lock(&_dist_mat_lock);
    omp_init_lock(&_hash_cache_lock);
    _n_threads = omp_get_max_threads();
    _hash_cache = CountingHashCache(_n_threads*2 + 1);
    // gets = 0;
    // loads = 0;
}

DistanceCalc::
~DistanceCalc()
{
    omp_set_lock(&_dist_mat_lock);
    if (_dist_mat == NULL) {
        omp_unset_lock(&_dist_mat_lock);
        return;
    }
    for (size_t i = 0; i < _n_samples; i++) {
        delete[] _dist_mat[i];
    }
    delete[] _dist_mat;
    omp_unset_lock(&_dist_mat_lock);
    omp_destroy_lock(&_dist_mat_lock);
}

float
DistanceCalc::
distance(khmer::CountingHash &a, khmer::CountingHash &b)
{
    return 0.0;
}


void
DistanceCalc::
calculate_pairwise(std::vector<std::string> &hash_fnames)
{
    // Create the distance matrix
    omp_set_lock(&_dist_mat_lock);
    _n_samples = hash_fnames.size();
    _dist_mat = new float *[_n_samples];
    for (size_t i = 0; i < _n_samples; i++) {
        _dist_mat[i] = new float[_n_samples];
    }
    omp_unset_lock(&_dist_mat_lock);

    #pragma omp parallel for schedule(dynamic) num_threads(_n_threads)
    for (size_t i = 0; i < _n_samples; i++) {
        CountingHashShrPtr ht1 = _get_hash(hash_fnames[i]);
        for (size_t j = 0; j < _n_samples; j++) {
            float dist = 0.0;
            if (i > j) {
                // Skip calculating the bottom half of the matrix
                continue;
            }
            CountingHashShrPtr ht2 = _get_hash(hash_fnames[j]);
            dist = this->distance(*ht1, *ht2);
            // Fill in both halves of the matrix
            _dist_mat[i][j] = dist;
            _dist_mat[j][i] = dist;
            #pragma omp critical
            {
                std::cerr << i << " x " << j << " done!" << std::endl;
            }
        }
    }
    std::cerr << "Done all!" << std::endl;
    //std::cerr << "Got " << gets - loads << " from cache, loaded " << loads << std::endl;
}

void
DistanceCalc::
set_sample_names(std::vector<std::string> &sample_names)
{
    for (auto name: sample_names) {
        _sample_names.push_back(name);
    }
}

void
DistanceCalc::
set_num_threads(int n_threads)
{
    _n_threads = n_threads;
}

void
DistanceCalc::
print_dist_mat(std::ostream &outstream)
{
    if (_dist_mat == NULL) {
        throw std::runtime_error("No distance matrix exists");
    }
    if (_sample_names.empty()) {
        for (size_t i = 0; i < _n_samples; i++) {
            _sample_names.push_back(std::to_string(i));
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
            outstream << "\t" << _dist_mat[i][j];
        }
        outstream << std::endl;
    }
}

CountingHashShrPtr
DistanceCalc::
_get_hash(std::string &filename)
{
    CountingHashShrPtr ret;
    while (1) {
        try {
            ret = _hash_cache.get(filename);
            //__sync_fetch_and_add(&gets, 1);
            return ret;
        } catch (std::range_error &err) {
            omp_set_lock(&_hash_cache_lock);
            //loads++;
            CountingHashShrPtr ht = \
                    std::make_shared<khmer::CountingHash>(1, 1);
            ht->load(filename);
            _hash_cache.put(filename, ht);
            //std::cerr << "Loaded " << filename << std::endl;
            omp_unset_lock(&_hash_cache_lock);
        }
    }
}

void
DistanceCalc::
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
