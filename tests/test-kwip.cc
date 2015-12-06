/*
 * ============================================================================
 *
 *       Filename:  test-kwip.cc
 *    Description:  Tests of kwip
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */


#include "helpers.hh"
#include "kernels/wip.hh"

#include "catch.hpp"


TEST_CASE("Test kwip.kernel method", "[kernel]") {
    kwip::metrics::WIPKernel kernel;
    using CountingHash = khmer::CountingHash;
    std::ostringstream output;

    kernel.outstream = &output;

    SECTION("Calc pairwise distances between the defined set") {
        CountingHash a{1, 1};
        CountingHash b{1, 1};
        std::vector<std::string> filenames {
            "data/defined-1.ct",
            "data/defined-2.ct",
            "data/defined-3.ct",
        };

        kernel.calculate_pairwise(filenames);

        std::vector<std::vector<float>> kmat_expt {
            {15.611,      9.18296,   6.42807},
            {9.18296,     10.1013,   0.918296},
            {6.42807,     0.918296,  7.34637},
        };
        std::vector<std::vector<float>> dmat_expt {
            {0.0,       0.733113,   0.894153},
            {0.733113,  0.0,        1.33671},
            {0.894153,  1.33671,    0.0},
        };

        float **kmat = kernel.get_kernel_matrix();
        float **dmat = kernel.get_distance_matrix();

        for (size_t i = 0; i < kmat_expt.size(); i++) {
            for (size_t j = 0; j < kmat_expt.size(); j++) {
                CAPTURE(i);
                CAPTURE(j);
                CAPTURE(dmat[i][j]);
                CAPTURE(kmat[i][j]);

                CHECK(dmat[i][j] == Approx(dmat_expt[i][j]));
                CHECK(kmat[i][j] == Approx(kmat_expt[i][j]));
            }
        }
    }
}
