#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "HullToObj.h"
#include <Hull/Coordinate.h>

#include <algorithm>
#include <math.h>
#include <sstream>

namespace hull {
bool operator==(const hull::Coordinate &a, const hull::Coordinate &b) {
  return (a.x == b.x) && (a.y == b.y) && (a.z == b.z);
}
} // namespace hull

hull::Hull make_hull(const std::vector<hull::Coordinate> &vertices) {
  hull::Hull hull(vertices[0], vertices[1], vertices[2], vertices[3]);
  auto it_vertices = vertices.begin();
  std::advance(it_vertices, 4);
  std::for_each(
      it_vertices, vertices.end(),
      [&hull](const hull::Coordinate &vertex) { hull.update(vertex); });
  return hull;
}

bool check_normals(const hull::Hull &subject) {
  const auto &vertices = subject.getVertices();
  hull::Coordinate mid_point;
  mid_point.x =
      0.25f * (vertices[0].x + vertices[1].x + vertices[2].x + vertices[3].x);
  mid_point.y =
      0.25f * (vertices[0].y + vertices[1].y + vertices[2].y + vertices[3].y);
  mid_point.z =
      0.25f * (vertices[0].z + vertices[1].z + vertices[2].z + vertices[3].z);
  for (const auto &facet : subject.getFacets()) {
    hull::Coordinate delta;
    hull::diff(delta, vertices[facet.vertexA], mid_point);
    if (hull::dot(vertices[facet.vertexA], facet.normal) <= 0) {
      return false;
    }
  }
  return true;
}

TEST_CASE("Simple thetraedron") {
  std::vector<hull::Coordinate> vertices = {
      {0, 0, 0}, {5, 5, 0}, {5, -5, 0}, {5, 0, 5}};

  auto hull = make_hull(vertices);

  hull::toObj(hull, "Thetraedron.obj");

  CHECK(hull.getVertices() == vertices);
  CHECK(check_normals(hull));
}

// TEST_CASE("Invalid thetraedrons") {

// }

float to_rad(const float angle) { return angle * M_PI / 180.0; }

std::vector<float> linspace(const float min, const float max,
                            const std::size_t intervals) {
  const float delta = (max - min) / intervals;
  std::vector<float> result = {min};
  for (std::size_t k = 0; k < intervals; ++k) {
    result.push_back(result.back() + delta);
  }
  return result;
}

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
  std::size_t angular_samples = GENERATE(3);
  // std::size_t angular_samples = GENERATE(3, 5, 10, 20, 50, 100, 500);

  std::vector<hull::Coordinate> vertices = make_sphere_cloud(angular_samples);

  auto hull = make_hull(vertices);

  std::stringstream log_name;
  log_name << "SphereCloud-" << std::to_string(angular_samples) << ".obj";
  hull::toObj(hull, log_name.str());

  CHECK(hull.getVertices() == vertices);
  CHECK(check_normals(hull));
}
