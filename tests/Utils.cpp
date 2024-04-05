#include "Utils.h"
#include <algorithm>
#include <fstream>
#include <sstream>

namespace hull {
void toObj(const HullContext &hull, const std::filesystem::path &fileName) {
  std::ofstream stream(fileName);
  for (const auto &vertex : hull.vertices) {
    stream << 'v';
    stream << ' ' << vertex.x;
    stream << ' ' << vertex.y;
    stream << ' ' << vertex.z;
    stream << std::endl;
  }
  stream << std::endl;
  for (auto *facet : hull.faces) {
    stream << 'f';
    stream << ' ' << std::to_string(facet->vertexA + 1);
    stream << ' ' << std::to_string(facet->vertexB + 1);
    stream << ' ' << std::to_string(facet->vertexC + 1);
    stream << std::endl;
  }
}

Logger::Logger() {
  if (std::filesystem::exists(LOG_FOLDER)) {
    std::filesystem::remove_all(LOG_FOLDER);
  }
  std::filesystem::create_directories(LOG_FOLDER);
}

Logger &Logger::get() {
  static Logger res = Logger{};
  return res;
}

std::filesystem::path Logger::makeLogfileName(const std::string &topicName) {
  auto it = topics.find(topicName);
  if (it == topics.end()) {
    it = topics.emplace(topicName, 0).first;
    std::filesystem::create_directories(std::filesystem::path{LOG_FOLDER} /
                                        topicName);
  }
  std::size_t c = ++it->second;
  std::filesystem::path res{LOG_FOLDER};
  res /= topicName;
  res /= "log-" + std::to_string(c) + ".obj";
  return res;
}

bool check_normals(const HullContext &hull) {
  const auto &vertices = hull.vertices;
  hull::Coordinate mid_point;
  mid_point.x =
      0.25f * (vertices[0].x + vertices[1].x + vertices[2].x + vertices[3].x);
  mid_point.y =
      0.25f * (vertices[0].y + vertices[1].y + vertices[2].y + vertices[3].y);
  mid_point.z =
      0.25f * (vertices[0].z + vertices[1].z + vertices[2].z + vertices[3].z);
  for (const auto &facet : hull.faces) {
    hull::Coordinate delta;
    hull::diff(delta, vertices[facet->vertexA], mid_point);
    if (hull::dot(delta, facet->normal) <= 0) {
      return false;
    }
  }
  return true;
}

void StepsLogger::hullChanges(const Notification &notification) {
  check_updated_mesh(notification);
  auto path = Logger::get().makeLogfileName(log_name.string());
  hull::toObj(notification.context, path);
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
  // check normals
  if (!check_normals(notification.context)) {
    throw std::runtime_error{"Invalid normals after update"};
  }
  // check connectivity
  const auto &facets = notification.context.faces;
  const auto &vertices = notification.context.vertices;
  for (const auto *facet : facets) {
    if (facet == facet->neighbourAB || facet == facet->neighbourBC ||
        facet == facet->neighbourCA) {
      throw std::runtime_error{"Neighbour of facet pointing to itself"};
    }
    if (!is_connected(*facet->neighbourAB, facet) ||
        !is_connected(*facet->neighbourBC, facet) ||
        !is_connected(*facet->neighbourCA, facet)) {
      throw std::runtime_error{"Neighbour not connected to this facet"};
    }
  }
}
} // namespace hull

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
