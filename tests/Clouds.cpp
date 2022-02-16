#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "Utils.h"
#include <Hull/Coordinate.h>

#include <algorithm>
#include <math.h>

namespace hull {
bool are_same(const std::vector<hull::Coordinate> &a,
              const std::vector<hull::Coordinate> &b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (const auto &a_vertex : a) {
    if (std::find_if(b.begin(), b.end(),
                     [&a_vertex](const hull::Coordinate &b_vertex) {
                       hull::Coordinate diff;
                       hull::diff(diff, b_vertex, a_vertex);
                       return hull::norm(diff) < 1e-5;
                     }) == b.end()) {
    }
  }
  return true;
}
} // namespace hull

constexpr float GREEK_PI = 3.14159f;
float to_rad(const float angle) { return angle * GREEK_PI / 180.f; }

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

bool is_connected(const hull::Facet &subject, const std::size_t facet_to_find) {
  if (subject.neighbourAB == facet_to_find) {
    return true;
  }
  if (subject.neighbourBC == facet_to_find) {
    return true;
  }
  if (subject.neighbourCA == facet_to_find) {
    return true;
  }
  return false;
}

class StepsLogger : public hull::Observer {
public:
  StepsLogger(const std::string &log_name) : log_name(log_name){};

  void hullChanges(const Notification &notification) override {
    check_updated_mesh(notification);
    hull::toObj(notification.vertices, notification.facets,
                generate_obj_log_name(log_name));
  };

  void check_updated_mesh(const Notification &notification) const {
    // check connectivity
    for (std::size_t facet_id = 0; facet_id < notification.facets.size();
         ++facet_id) {
      if (facet_id == notification.facets[facet_id].neighbourAB) {
        throw std::runtime_error{"Neighbour of facet pointing to itself"};
      }
      if (!is_connected(
              notification.facets[notification.facets[facet_id].neighbourAB],
              facet_id)) {
        throw std::runtime_error{"Neighbour not connected to this facet"};
      }

      if (facet_id == notification.facets[facet_id].neighbourBC) {
        throw std::runtime_error{"Neighbour of facet pointing to itself"};
      }
      if (!is_connected(
              notification.facets[notification.facets[facet_id].neighbourBC],
              facet_id)) {
        throw std::runtime_error{"Neighbour not connected to this facet"};
      }

      if (facet_id == notification.facets[facet_id].neighbourCA) {
        throw std::runtime_error{"Neighbour of facet pointing to itself"};
      }
      if (!is_connected(
              notification.facets[notification.facets[facet_id].neighbourCA],
              facet_id)) {
        throw std::runtime_error{"Neighbour not connected to this facet"};
      }
    }
  }

private:
  const std::string log_name;
};

TEST_CASE("Sphere clouds") {
  std::size_t angular_samples = GENERATE(3);
  // std::size_t angular_samples = GENERATE(3, 5, 10, 20, 50, 100, 500);

  std::vector<hull::Coordinate> vertices = make_sphere_cloud(angular_samples);

  hull::Hull hull(vertices[0], vertices[1], vertices[vertices.size() - 2],
                  vertices.back());

  StepsLogger logger(std::string("SphereCloud-") +
                     std::to_string(angular_samples));
  hull.setObserver(logger);

  auto it_vertices = vertices.begin();
  std::advance(it_vertices, 4);
  std::for_each(
      it_vertices, vertices.end(),
      [&hull](const hull::Coordinate &vertex) { hull.update(vertex); });

  CHECK(are_same(hull.getVertices(), vertices));
  CHECK(check_normals(hull));
}
