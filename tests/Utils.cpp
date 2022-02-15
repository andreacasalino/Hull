#include "Utils.h"
#include <fstream>
#include <sstream>

namespace hull {
void toObj(const Hull &subject, const std::string &fileName) {
  std::ofstream stream(fileName);
  for (const auto &vertex : subject.getVertices()) {
    stream << 'v';
    stream << ' ' << vertex.x;
    stream << ' ' << vertex.y;
    stream << ' ' << vertex.z;
    stream << std::endl;
  }
  stream << std::endl;
  for (const auto &facet : subject.getFacets()) {
    stream << 'f';
    stream << ' ' << std::to_string(facet.vertexA + 1);
    stream << ' ' << std::to_string(facet.vertexB + 1);
    stream << ' ' << std::to_string(facet.vertexC + 1);
    stream << std::endl;
  }
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
    if (hull::dot(delta, facet.normal) <= 0) {
      return false;
    }
  }
  return true;
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
