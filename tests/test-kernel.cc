/*
 * ============================================================================
 *
 *       Filename:  test-kernel.cc
 *    Description:  Tests of the core 'kernel.cc' functionality
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include "catch.hpp"
#include "helpers.hh"

#include "kernel.hh"

TEST_CASE("Test kernel before computation", "[kernel]") {
    kmerclust::Kernel kernel;

    SECTION("Test num_threads") {
        // Sensible default
        REQUIRE(kernel.num_threads == omp_get_max_threads());

        // Explicit setting
        kernel.num_threads = 5;
        REQUIRE(kernel.num_threads == 5);
    }

    SECTION("Test verbosity") {
        // Sensible default
        REQUIRE(kernel.verbosity == 1);

        // Explicit setting
        kernel.verbosity = 2;
        REQUIRE(kernel.verbosity == 2);
    }

    SECTION("Test name") {
        REQUIRE(kernel.name == "Base Class");
    }

    SECTION("Test blurb") {
        // Don't test content
        REQUIRE(kernel.blurb.size() > 0);
    }

    SECTION("Test sample_names") {
        REQUIRE(kernel.sample_names.empty());
        REQUIRE(kernel.num_samples == 0);
    }

    SECTION("Test print_kernel_mat") {
        REQUIRE_THROWS_AS(kernel.print_kernel_mat(), std::runtime_error);
    }

    SECTION("Test print_dist_mat") {
        REQUIRE_THROWS_AS(kernel.print_kernel_mat(), std::runtime_error);
    }
}


TEST_CASE("Test kernel.kernel method", "[kernel]") {
    using CountingHash = khmer::CountingHash;
    kmerclust::Kernel kernel;

    SECTION("Same hash dimensions") {
        CountingHash a{10, 10000};
        CountingHash b{10, 10000};

        REQUIRE(kernel.kernel(a, b) == 0.0);
    }

    SECTION("Different table size") {
        CountingHash a{10, 10000};
        CountingHash b{10, 1000};

        REQUIRE_THROWS_AS(kernel.kernel(a, b), std::runtime_error);
    }

    SECTION("Different K") {
        CountingHash a{10, 10000};
        CountingHash b{12, 10000};

        REQUIRE_THROWS_AS(kernel.kernel(a, b), std::runtime_error);
    }

    SECTION("Different K and table size") {
        CountingHash a{10, 10000};
        CountingHash b{12, 1000};

        REQUIRE_THROWS_AS(kernel.kernel(a, b), std::runtime_error);
    }
}

