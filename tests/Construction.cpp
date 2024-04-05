#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "Utils.h"
#include <Hull/Coordinate.h>

TEST_CASE("Hull construction") {
  SECTION("expected to not throw tests") {
    auto vertices = GENERATE(
        std::vector<hull::Coordinate>{
            {0, 0, 0}, {5, 5, 0}, {5, -5, 0}, {5, 0, 5}},
        std::vector<hull::Coordinate>{
            {0, 0, -1}, {0, 0, 1}, {1, 1, 0}, {1, -1, 0}});

    hull::Hull hull(vertices[0], vertices[1], vertices[2], vertices[3]);

    hull::toObj(hull.getContext(),
                hull::Logger::get().makeLogfileName("HullConstruction"));

    CHECK(hull::check_normals(hull.getContext()));
  }

  SECTION("expected to throw tests") {
    auto vertices = GENERATE(
        std::vector<hull::Coordinate>{
            {0, 0, 0}, {5, 5, 0}, {5, -5, 0}, {5, 0, 0}},
        std::vector<hull::Coordinate>{
            {-1, 0, 0}, {0, 0, 0}, {1, 1, 0}, {1, -1, 0}});

    CHECK_THROWS(
        hull::Hull{vertices[0], vertices[1], vertices[2], vertices[3]});
  }
}