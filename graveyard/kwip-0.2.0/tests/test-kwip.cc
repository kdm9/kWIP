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

#include <Eigen/Core>
using Eigen::Matrix3d;
using Eigen::MatrixXd;

#include "helpers.hh"
#include "kernels/wip.hh"

#include "catch.hpp"


TEST_CASE("Test kwip.kernel method", "[kernel]") {
    kwip::metrics::WIPKernel kernel;
    MatrixXd kmat, nkmat, dmat;
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

        MatrixXd kmat_expt(3, 3);
        kmat_expt <<
             15.611,      9.18296,   6.42807,
             9.18296,     10.1013,   0.918296,
             6.42807,     0.918296,  7.34637;

        MatrixXd dmat_expt(3, 3);
        dmat_expt <<
             0.0,       0.733113,   0.894153,
             0.733113,  0.0,        1.33671,
             0.894153,  1.33671,    0.0;

        kernel.get_kernel_matrix(kmat);
        kernel.get_norm_kernel_matrix(nkmat);
        kernel.get_distance_matrix(dmat);

        CHECK(kmat.isApprox(kmat_expt, 1e-4));
        CHECK(dmat.isApprox(dmat_expt, 1e-4));
    }
}
