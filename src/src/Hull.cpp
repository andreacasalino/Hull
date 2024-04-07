/**
 * Author:    Andrea Casalino
 * Created:   03.12.2019
 *
 * report any bug to andrecasa91@gmail.com.
 **/

#include <Hull/Definitions.h>
#include <Hull/Error.h>
#include <Hull/Hull.h>

#include <algorithm>
#include <array>
#include <stack>
#include <unordered_map>

namespace hull {
FacetAllocator::FacetAllocator() { buffers.emplace_back(10); }

FacetAllocator::Buffer::Buffer(std::size_t capacity)
    : capacity{capacity}, buffer{new char[sizeof(Facet) * capacity]} {}

FacetAllocator::Buffer::~Buffer() {
  Facet *ptr = reinterpret_cast<Facet *>(buffer);
  for (std::size_t k = 0; k < size; ++k, ++ptr) {
    ptr->~Facet(); // TODO given the current definition of Facet, is it actually
                   // needed?
  }
  delete[] buffer;
}

Facet *FacetAllocator::getFacet() {
  if (!outstanding.empty()) {
    Facet *res = outstanding.front();
    outstanding.pop_front();
    return res;
  }

  Buffer *last_buffer = &buffers.back();
  auto create = [&last_buffer]() {
    Facet *ptr = reinterpret_cast<Facet *>(last_buffer->buffer);
    ptr += last_buffer->size;
    new (ptr) Facet{};
    ++last_buffer->size;
    return ptr;
  };

  if (last_buffer->size < last_buffer->capacity) {
    return create();
  }
  last_buffer = &buffers.emplace_back(last_buffer->capacity * 2);
  return create();
}

Hull::Hull(const Coordinate &A, const Coordinate &B, const Coordinate &C,
           const Coordinate &D) {
  this->initThetraedron(A, B, C, D);
};

Hull::Hull(const Coordinate &A, const Coordinate &B, const Coordinate &C,
           const Coordinate &D, Observer &obs)
    : observer(&obs) {
  this->initThetraedron(A, B, C, D);
}

void Hull::setObserver(Observer &obs) { this->observer = &obs; }

namespace {
void recomputeNormal(Facet &subject, const HullContext &ctxt) {
  Coordinate delta1, delta2;

  const auto &getPoint = [&vertices = ctxt.vertices](std::size_t index) {
    return vertices[index];
  };

  diff(delta1, getPoint(subject.vertexA), getPoint(subject.vertexC));
  diff(delta2, getPoint(subject.vertexB), getPoint(subject.vertexC));

  cross(subject.normal, delta1, delta2);

  diff(delta1, ctxt.internalPoint, getPoint(subject.vertexA));
  float dot_normal = dot(subject.normal, delta1);
  if (dot_normal >= 0.f) {
    invert(subject.normal);
  }
  // normalize
  dot_normal = norm(subject.normal);
  if (dot_normal < 1e-7) {
    prod(subject.normal, COEFF_NORMAL_DIRECTION);
  } else {
    dot_normal = 1.f / dot_normal;
    prod(subject.normal, dot_normal);
  }
}
} // namespace

Facet *Hull::makeFacet(std::size_t vertexA, std::size_t vertexB,
                       std::size_t vertexC) {
  Facet *result = allocator.getFacet();
  result->vertexA = vertexA;
  result->vertexB = vertexB;
  result->vertexC = vertexC;
  recomputeNormal(*result, context);
  return result;
}

void Hull::initThetraedron(const Coordinate &A, const Coordinate &B,
                           const Coordinate &C, const Coordinate &D) {
  {
    // check the thetraedron has a non zero volume
    Coordinate delta_1;
    diff(delta_1, D, A);
    Coordinate delta_2;
    diff(delta_2, D, B);
    Coordinate delta_3;
    diff(delta_3, D, C);

    Coordinate cross_1_2;
    cross(cross_1_2, delta_1, delta_2);
    if (abs(dot(cross_1_2, delta_3)) < HULL_GEOMETRIC_TOLLERANCE) {
      throw Error("intial thetraedron volume too small: make sure the convex "
                  "shape is actually a 3d shape");
    }
  }

  // computation of the midpoint of the thetraedron
  this->context.internalPoint.x = 0.25f * (A.x + B.x + C.x + D.x);
  this->context.internalPoint.y = 0.25f * (A.y + B.y + C.y + D.y);
  this->context.internalPoint.z = 0.25f * (A.z + B.z + C.z + D.z);

  // build the tethraedron
  // ABC->0; ABD->1; ACD->2; BCD->3
  auto &vertices = context.vertices;
  vertices.push_back(A);
  vertices.push_back(B);
  vertices.push_back(C);
  vertices.push_back(D);
  std::array<Facet *, 4> facets{makeFacet(0, 1, 2), makeFacet(0, 1, 3),
                                makeFacet(0, 2, 3), makeFacet(1, 2, 3)};
  // setup initial connectivity
  // ABC
  facets[0]->neighbourAB = facets[1];
  facets[0]->neighbourBC = facets[3];
  facets[0]->neighbourCA = facets[2];
  // ABD
  facets[1]->neighbourAB = facets[0];
  facets[1]->neighbourBC = facets[3];
  facets[1]->neighbourCA = facets[2];
  // ACD
  facets[2]->neighbourAB = facets[0];
  facets[2]->neighbourBC = facets[3];
  facets[2]->neighbourCA = facets[1];
  // BCD
  facets[3]->neighbourAB = facets[0];
  facets[3]->neighbourBC = facets[2];
  facets[3]->neighbourCA = facets[1];
  context.faces.insert(facets.begin(), facets.end());

  if (this->observer) {
    Observer::Notification update_notif{context};
    update_notif.added.insert(update_notif.added.end(), facets.begin(),
                              facets.end());
    observer->hullChanges(update_notif);
  }
}

namespace {
float facet_point_distance(const Facet &facet, const Coordinate &point,
                           const HullContext &ctxt) {
  float distance = facet.normal.x * (point.x - ctxt.vertices[facet.vertexA].x);
  distance += facet.normal.y * (point.y - ctxt.vertices[facet.vertexA].y);
  distance += facet.normal.z * (point.z - ctxt.vertices[facet.vertexA].z);
  return distance;
}
} // namespace

void Hull::update(const Coordinate &vertex_of_new_cone) {
  // find first visible facet
  auto it = std::find_if(
      context.faces.begin(), context.faces.end(), [&](Facet *facet) {
        return HULL_GEOMETRIC_TOLLERANCE <
               facet_point_distance(*facet, vertex_of_new_cone, context);
      });
  if (it == context.faces.end()) {
    throw Error{"Vertex passed to update hull is not outside it"};
  }
  this->update_(vertex_of_new_cone, *it);
}

void Hull::update(const Coordinate &vertex_of_new_cone,
                  Facet *starting_facet_for_expansion) {
  const auto &faces = context.faces;
  if (HULL_GEOMETRIC_TOLLERANCE <=
          facet_point_distance(*starting_facet_for_expansion,
                               vertex_of_new_cone, context) ||
      context.faces.find(starting_facet_for_expansion) == context.faces.end()) {
    throw Error{"The specified starting facet is not valid"};
  }
  this->update_(vertex_of_new_cone, starting_facet_for_expansion);
}

namespace {
enum ConnectivityCase { AB, BC, CA };

ConnectivityCase find_connectivity_case(const Facet &subject,
                                        const Facet *neighbour_to_find) {
  if (neighbour_to_find == subject.neighbourAB) {
    return ConnectivityCase::AB;
  }
  if (neighbour_to_find == subject.neighbourBC) {
    return ConnectivityCase::BC;
  }
  if (neighbour_to_find == subject.neighbourCA) {
    return ConnectivityCase::CA;
  }
  throw Error{"neighbour index not found"};
}

struct Edges {
  using Info = std::pair<Facet *, ConnectivityCase>;

  void add(std::size_t vertexFirst, std::size_t vertexSecond,
           const Info &to_add) {
    std::size_t index = edges.size();
    edges.emplace_back(std::make_tuple(vertexFirst, vertexSecond, to_add));
    edges_table[vertexFirst][vertexSecond] = std::make_pair(index, to_add);
    edges_table[vertexSecond][vertexFirst] = std::make_pair(index, to_add);
  }

  // return the other edge having vertexA as vertex, but no vertexB as the
  // second one
  std::size_t connectedEdgeIndex(std::size_t vertexFirst,
                                 std::size_t vertexSecond) const {
    auto it_primary = edges_table.find(vertexFirst);
    auto it_one = it_primary->second.begin();
    auto it_second = it_one;
    ++it_second;
    return it_one->first == vertexSecond ? it_second->second.first
                                         : it_one->second.first;
  }

  std::vector<std::tuple<std::size_t, std::size_t, Info>> edges;

private:
  std::unordered_map<
      std::size_t,
      std::unordered_map<std::size_t, std::pair<std::size_t, Info>>>
      edges_table;
};

struct VisibleCone {
  Edges edges;
  std::unordered_set<Facet *> visible_faces;
};

VisibleCone computeVisibleCone(const Coordinate &vertex_of_new_cone,
                               Facet *starting_facet, const HullContext &ctxt) {
  auto &facets = ctxt.faces;
  const auto &vertices = ctxt.vertices;

  VisibleCone res;
  auto &edges = res.edges;
  auto &visible_group = res.visible_faces;

  // this set contains facets we know that are visible, but
  // whose neighbouring was not already visited
  std::stack<Facet *> facets_to_visit;
  facets_to_visit.push(starting_facet);
  while (!facets_to_visit.empty()) {
    auto *to_visit = facets_to_visit.top();
    facets_to_visit.pop();
    if (visible_group.find(to_visit) != visible_group.end()) {
      // already visited ...
      continue;
    }
    visible_group.emplace(to_visit);

    auto process_neighbour = [&](Facet *neighbour, std::size_t vertexFirst,
                                 std::size_t vertexSecond) {
      if (facet_point_distance(*neighbour, vertex_of_new_cone, ctxt) >
          HULL_GEOMETRIC_TOLLERANCE) {
        facets_to_visit.push(neighbour);
      } else {
        Edges::Info info = std::make_pair(
            neighbour, find_connectivity_case(*neighbour, to_visit));
        edges.add(vertexFirst, vertexSecond, info);
      }
    };

    process_neighbour(to_visit->neighbourAB, to_visit->vertexA,
                      to_visit->vertexB);
    process_neighbour(to_visit->neighbourBC, to_visit->vertexB,
                      to_visit->vertexC);
    process_neighbour(to_visit->neighbourCA, to_visit->vertexC,
                      to_visit->vertexA);
  }
  return res;
}
} // namespace

void Hull::update_(const Coordinate &vertex_of_new_cone,
                   Facet *starting_facet_for_expansion) {
  auto visibility_cone = computeVisibleCone(
      vertex_of_new_cone, starting_facet_for_expansion, context);
  Observer::Notification update_notif{context};

  std::vector<Facet *> facets_cone;
  std::size_t cone_size = visibility_cone.edges.edges.size();
  if (visibility_cone.visible_faces.size() <= cone_size) {
    // more edges than visible facets
    facets_cone.insert(facets_cone.end(), visibility_cone.visible_faces.begin(),
                       visibility_cone.visible_faces.end());
    while (facets_cone.size() < cone_size) {
      auto *f = allocator.getFacet();
      facets_cone.push_back(f);
      context.faces.insert(f);
    }
    if (observer) {
      update_notif.changed.insert(update_notif.changed.end(),
                                  visibility_cone.visible_faces.begin(),
                                  visibility_cone.visible_faces.end());
      update_notif.added.insert(update_notif.added.end(),
                                facets_cone.begin() +
                                    visibility_cone.visible_faces.size(),
                                facets_cone.end());
    }
  } else {
    // more visible facets than edges
    auto it = visibility_cone.visible_faces.begin();
    for (std::size_t c = 0; c < cone_size; ++c, ++it) {
      facets_cone.push_back(*it);
    }
    for (std::size_t c = cone_size; c < visibility_cone.visible_faces.size();
         ++c, ++it) {
      allocator.markAsRecyclable(*it);
      context.faces.erase(*it);
    }
    if (observer) {
      it = visibility_cone.visible_faces.begin();
      for (std::size_t c = 0; c < cone_size; ++c, ++it) {
        update_notif.changed.push_back(*it);
      }
      for (std::size_t c = cone_size; c < visibility_cone.visible_faces.size();
           ++c, ++it) {
        update_notif.removed.push_back(*it);
      }
    }
  }
  std::size_t added_vertex = context.vertices.size();
  context.vertices.emplace_back(vertex_of_new_cone);

  // build cone of new facets
  const auto &edges = visibility_cone.edges.edges;
  for (std::size_t k = 0; k < cone_size; ++k) {
    const auto &[vertexFirst, vertexSecond, info] = edges[k];
    auto *facet_to_build = facets_cone[k];
    facet_to_build->vertexA = vertexFirst;
    facet_to_build->vertexB = vertexSecond;
    facet_to_build->vertexC = added_vertex;

    facet_to_build->neighbourAB = info.first;
    switch (info.second) {
    case ConnectivityCase::AB:
      facet_to_build->neighbourAB->neighbourAB = facet_to_build;
      break;
    case ConnectivityCase::BC:
      facet_to_build->neighbourAB->neighbourBC = facet_to_build;
      break;
    case ConnectivityCase::CA:
      facet_to_build->neighbourAB->neighbourCA = facet_to_build;
      break;
    }

    facet_to_build->neighbourCA =
        facets_cone[visibility_cone.edges.connectedEdgeIndex(vertexFirst,
                                                             vertexSecond)];

    facet_to_build->neighbourBC =
        facets_cone[visibility_cone.edges.connectedEdgeIndex(vertexSecond,
                                                             vertexFirst)];

    recomputeNormal(*facet_to_build, context);
  }

  // update mid point in order to have more precise future iterations
  {
    const auto &vertices = context.vertices;
    prod(context.internalPoint, static_cast<float>(vertices.size() - 1));
    context.internalPoint.x += vertices.back().x;
    context.internalPoint.y += vertices.back().y;
    context.internalPoint.z += vertices.back().z;
    prod(context.internalPoint, 1.f / static_cast<float>(vertices.size()));
  }

  if (this->observer) {
    observer->hullChanges(update_notif);
  }
}
} // namespace hull
