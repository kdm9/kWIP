/*
 * ============================================================================
 *
 *       Filename:  test-lrucache.cc
 *    Description:  Test of the LRUcache implementation
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include "catch.hpp"
#include "helpers.hh"

#include "lrucache.hpp"

TEST_CASE("LRU Cache basic operations", "[lru-cache]") {
    cache::lru_cache<int, int> lru(3);

    SECTION("cache.put overwrites current value") {
        lru.put(1, 1);
        REQUIRE_NOTHROW(lru.put(1, 4));
        REQUIRE(lru.get(1) == 4);
    }

    SECTION("cache.get retrives correct values") {
        lru.put(1, 1);
        REQUIRE(lru.get(1) == 1);
        REQUIRE_THROWS_AS(lru.get(5), std::range_error);
    }

    SECTION("cache.put expunges least-recently added element if unused") {
        lru.put(1, 1);
        lru.put(2, 2);
        lru.put(3, 3);
        lru.put(4, 4);
        REQUIRE(lru.get(2) == 2);
        REQUIRE(lru.get(3) == 3);
        REQUIRE(lru.get(4) == 4);
        REQUIRE_THROWS_AS(lru.get(1), std::range_error);
    }

    SECTION("cache.put expunges least-recently used element") {
        lru.put(1, 1);
        lru.put(2, 2);
        lru.put(3, 3);
        REQUIRE(lru.get(2) == 2);
        REQUIRE(lru.get(3) == 3);
        REQUIRE(lru.get(1) == 1);
        // The least recently used element is 2
        lru.put(4, 4);
        // Test 2, and only 2, is expunged
        REQUIRE_THROWS_AS(lru.get(2), std::range_error);
        REQUIRE(lru.get(1) == 1);
        REQUIRE(lru.get(3) == 3);
        REQUIRE(lru.get(4) == 4);
    }

    SECTION("cache.exists is accurate") {
        lru.put(1, 1);
        REQUIRE(lru.exists(1));
        REQUIRE(!lru.exists(2));
    }

    SECTION("cache.size is accurate") {
        for (int i = 1; i < 4; i++) {
            lru.put(i, i);
            REQUIRE(lru.size() == i);
        }
    }
}
