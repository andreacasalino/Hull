#include "Utils.h"
#include <algorithm>
#include <fstream>
#include <sstream>

namespace hull {
void toObj(const std::vector<Coordinate> &vertices, const Facets &facets,
           const std::string &fileName) {
  std::ofstream stream(fileName);
  for (const auto &vertex : vertices) {
    stream << 'v';
    stream << ' ' << vertex.x;
    stream << ' ' << vertex.y;
    stream << ' ' << vertex.z;
    stream << std::endl;
  }
  stream << std::endl;
  for (const auto &facet : facets) {
    stream << 'f';
    stream << ' ' << std::to_string(facet->vertexA + 1);
    stream << ' ' << std::to_string(facet->vertexB + 1);
    stream << ' ' << std::to_string(facet->vertexC + 1);
    stream << std::endl;
  }
}

bool check_normals(const std::vector<Coordinate> &vertices,
                   const Facets &facets) {
  hull::Coordinate mid_point;
  mid_point.x =
      0.25f * (vertices[0].x + vertices[1].x + vertices[2].x + vertices[3].x);
  mid_point.y =
      0.25f * (vertices[0].y + vertices[1].y + vertices[2].y + vertices[3].y);
  mid_point.z =
      0.25f * (vertices[0].z + vertices[1].z + vertices[2].z + vertices[3].z);
  for (const auto &facet : facets) {
    hull::Coordinate delta;
    hull::diff(delta, vertices[facet->vertexA], mid_point);
    if (hull::dot(delta, facet->normal) <= 0) {
      return false;
    }
  }
  return true;
}

void StepsLogger::hullChanges(Notification &&notification) {
  check_updated_mesh(notification);
  hull::toObj(notification.context.vertices, notification.context.faces,
              generate_obj_log_name(log_name));
};

namespace {
bool is_connected(const hull::Facet &subject, const Facet *facet_to_find) {
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
} // namespace

void StepsLogger::check_updated_mesh(const Notification &notification) const {
  const auto &facets = notification.context.faces;
  const auto &vertices = notification.context.vertices;
  // check normals
  if (!check_normals(vertices, facets)) {
    throw std::runtime_error{"Invalid normals after update"};
  }
  // check connectivity
  for (const auto &facet : facets) {
    if (facet.get() == facet->neighbourAB) {
      throw std::runtime_error{"Neighbour of facet pointing to itself"};
    }
    if (!is_connected(*facet->neighbourAB, facet.get())) {
      throw std::runtime_error{"Neighbour not connected to this facet"};
    }

    if (facet.get() == facet->neighbourBC) {
      throw std::runtime_error{"Neighbour of facet pointing to itself"};
    }
    if (!is_connected(*facet->neighbourBC, facet.get())) {
      throw std::runtime_error{"Neighbour not connected to this facet"};
    }

    if (facet.get() == facet->neighbourCA) {
      throw std::runtime_error{"Neighbour of facet pointing to itself"};
    }
    if (!is_connected(*facet->neighbourCA, facet.get())) {
      throw std::runtime_error{"Neighbour not connected to this facet"};
    }
  }
}
} // namespace hull

#include <mutex>
#include <unordered_map>

namespace {
static std::mutex counters_mtx;
static std::unordered_map<std::string, std::size_t> counters;
} // namespace

std::string generate_obj_log_name(const std::string &file_name) {
  std::scoped_lock lock(counters_mtx);
  auto counters_it = counters.find(file_name);
  if (counters_it == counters.end()) {
    counters_it = counters.emplace(file_name, 0).first;
  } else {
    ++counters_it->second;
  }
  std::stringstream stream;
  stream << counters_it->second << '-' << file_name << ".obj";
  return stream.str();
}

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

hull::Hull fill_hull(const std::vector<hull::Coordinate> &vertices,
                     hull::Observer *obs) {
  if (vertices.size() < 4) {
    throw std::runtime_error{"invalid set of vertices"};
  }

  hull::Hull result(vertices.front(), vertices[1],
                    vertices[vertices.size() - 2], vertices.back());

  if (nullptr != obs) {
    result.setObserver(*obs);
  }

  auto vertices_begin = vertices.begin();
  std::advance(vertices_begin, 2);
  auto vertices_end = vertices.end();
  std::advance(vertices_end, -2);
  std::for_each(
      vertices_begin, vertices_end,
      [&result](const hull::Coordinate &vertex) { result.update(vertex); });

  return result;
}
