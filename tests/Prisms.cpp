#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "Utils.h"
#include <Hull/Coordinate.h>

std::vector<hull::Coordinate> make_prism(const std::size_t polygon_size,
                                         const float height) {
  std::vector<float> angles = linspace(0, to_rad(360), polygon_size);
  angles.pop_back();

  std::vector<hull::Coordinate> result;
  result.reserve(angles.size() * 2);
  for (const auto &angle : angles) {
    result.push_back(
        hull::Coordinate{10.f * cos(angle), 10.f * sin(angle), -height});
  }
  for (const auto &angle : angles) {
    result.push_back(
        hull::Coordinate{10.f * cos(angle), 10.f * sin(angle), height});
  }
  return result;
}

TEST_CASE("Prisms") {
  SECTION("normal height") {
    std::size_t polygon_size = GENERATE(3, 5, 10);

    const float h = 5;
    std::vector<hull::Coordinate> vertices = make_prism(polygon_size, h);

    hull::StepsLogger logger(std::string("Prism-h") + std::to_string(h) +
                             std::string("-") + std::to_string(polygon_size));

    auto hull = fill_hull(vertices, &logger);

    CHECK(are_same(hull.getVertices(), vertices));
  }

  SECTION("small height") {
    std::size_t polygon_size = GENERATE(3, 5, 10);

    const float h = 1;
    std::vector<hull::Coordinate> vertices = make_prism(polygon_size, h);

    hull::StepsLogger logger(std::string("Prism-h") + std::to_string(h) +
                             std::string("-") + std::to_string(polygon_size));

    auto hull = fill_hull(vertices, &logger);

    CHECK(are_same(hull.getVertices(), vertices));
  }
}
