#include "HullToObj.h"
#include <fstream>

namespace hull {
void toObj(const Hull &subject, const std::string &fileName) {
  std::ofstream stream(fileName);
  for (const auto &vertex : subject.getVertices()) {
    stream << 'v';
    stream << ' ' << vertex.x;
    stream << ' ' << vertex.y;
    stream << ' ' << vertex.z;
  }
  stream << std::endl;
  for (const auto &facet : subject.getFacets()) {
    stream << 'f';
    stream << ' ' << std::to_string(facet.vertexA + 1);
    stream << ' ' << std::to_string(facet.vertexB + 1);
    stream << ' ' << std::to_string(facet.vertexC + 1);
  }
}
} // namespace hull
