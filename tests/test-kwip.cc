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
#include "kernels/d2ent.hh"

#include "catch.hpp"


TEST_CASE("Test kwip.kernel method", "[kernel]") {
    kmerclust::metrics::KernelD2Ent kernel;
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
            {0.0286967307329178, 0.0286967437714338, 0.0286967437714338},
            {0.0286967437714338, 0.0573934875428677, 0.0},
            {0.0286967437714338, 0.0,                0.0573934875428677}
        };
        std::vector<std::vector<float>> dmat_expt {
            {0.0,               0.765366673469543, 0.765366673469543},
            {0.765366673469543, 0.0,               1.41421353816986},
            {0.765366673469543, 1.41421353816986,  0.0}
        };

        float **kmat = kernel.get_kernel_matrix();
        float **dmat = kernel.get_distance_matrix();

        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
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
