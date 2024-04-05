#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "Utils.h"
#include <Hull/Coordinate.h>

std::vector<hull::Coordinate>
make_sphere_cloud(const std::size_t angular_samples) {
  std::vector<float> angles1 =
      linspace(to_rad(-85), to_rad(85), angular_samples);
  std::vector<float> angles2 =
      linspace(to_rad(0), to_rad(350), angular_samples);

  std::vector<hull::Coordinate> result;
  result.reserve(angles1.size() * angles2.size());
  for (const auto &angle1 : angles1) {
    for (const auto &angle2 : angles2) {
      result.push_back(hull::Coordinate{10.f * cos(angle1) * cos(angle2),
                                        10.f * sin(angle1) * cos(angle2),
                                        10.f * sin(angle2)});
    }
  }
  return result;
}

TEST_CASE("Sphere clouds") {
  std::size_t angular_samples = 10;
  // std::size_t angular_samples = GENERATE(3, 5, 10, 20, 50);

  std::vector<hull::Coordinate> vertices = make_sphere_cloud(angular_samples);

  hull::StepsLogger logger(std::string("SphereCloud-") +
                           std::to_string(angular_samples));

  auto hull = fill_hull(vertices, &logger);

  CHECK(same_vertices(hull.getContext().vertices, vertices));
}
