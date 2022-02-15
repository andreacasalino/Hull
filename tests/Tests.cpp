#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "HullToObj.h"
#include <Hull/Coordinate.h>

#include <algorithm>
#include <sstream>

namespace hull {
bool operator==(const hull::Coordinate &a, const hull::Coordinate &b) {
  return (a.x == b.x) && (a.y == b.y) && (a.z == b.z);
}
} // namespace hull

hull::Hull make_hull(const std::vector<hull::Coordinate> &vertices) {
  hull::Hull hull(vertices[0], vertices[1], vertices[12], vertices[3]);
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

// TEST_CASE("Sphere clouds") {
//   std::size_t angles_samples = GENERATE(5, 10, 20, 50, 100, 500);

//   std::vector<hull::Coordinate> vertices = make_sphere_cloud(samples);

//   auto hull = make_hull(vertices);

//   hull::toObj(hull, "Thetraedron.obj");

//   std::stringstream log_name;
//   log_name << "SphereCloud-" << std::to_string(samples) << ".obj";
//   hull::toObj(hull, log_name.str());

//   CHECK(hull.getVertices() == vertices);
//   CHECK(check_normals(hull));
// }
