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
#include <limits>
#include <list>
#include <set>

namespace hull {
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

void Hull::recomputeNormal(Facet &subject) const {
  Coordinate delta1, delta2;
  const auto &vertices = vertices_and_faces.vertices;
  diff(delta1, vertices[subject.vertexA], vertices[subject.vertexC]);
  diff(delta2, vertices[subject.vertexB], vertices[subject.vertexC]);

  cross(subject.normal, delta1, delta2);

  diff(delta1, this->Mid_point, vertices[subject.vertexA]);
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

FacetPtr Hull::makeFacet(const std::size_t vertexA, const std::size_t vertexB,
                         const std::size_t vertexC) const {
  FacetPtr result = std::make_unique<Facet>();
  result->vertexA = vertexA;
  result->vertexB = vertexB;
  result->vertexC = vertexC;
  recomputeNormal(*result);
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
  this->Mid_point.x = 0.25f * (A.x + B.x + C.x + D.x);
  this->Mid_point.y = 0.25f * (A.y + B.y + C.y + D.y);
  this->Mid_point.z = 0.25f * (A.z + B.z + C.z + D.z);

  // build the tethraedron
  // ABC->0; ABD->1; ACD->2; BCD->3
  auto &vertices = vertices_and_faces.vertices;
  auto &facets = vertices_and_faces.faces;
  vertices.reserve(4);
  vertices.push_back(A);
  vertices.push_back(B);
  vertices.push_back(C);
  vertices.push_back(D);
  facets.reserve(4);
  facets.emplace_back(makeFacet(0, 1, 2));
  facets.emplace_back(makeFacet(0, 1, 3));
  facets.emplace_back(makeFacet(0, 2, 3));
  facets.emplace_back(makeFacet(1, 2, 3));
  // setup initial connectivity
  // ABC
  facets[0]->neighbourAB = facets[1].get();
  facets[0]->neighbourBC = facets[3].get();
  facets[0]->neighbourCA = facets[2].get();
  // ABD
  facets[1]->neighbourAB = facets[0].get();
  facets[1]->neighbourBC = facets[3].get();
  facets[1]->neighbourCA = facets[2].get();
  // ACD
  facets[2]->neighbourAB = facets[0].get();
  facets[2]->neighbourBC = facets[3].get();
  facets[2]->neighbourCA = facets[1].get();
  // BCD
  facets[3]->neighbourAB = facets[0].get();
  facets[3]->neighbourBC = facets[2].get();
  facets[3]->neighbourCA = facets[1].get();

  if (nullptr != this->observer) {
    observer->hullChanges(Observer::Notification{
        {},
        {facets[0].get(), facets[1].get(), facets[2].get(), facets[3].get()},
        {},
        vertices_and_faces});
  }
}

namespace {
float facet_point_distance(const std::vector<Coordinate> &vertices,
                           const Facet &facet, const Coordinate &point) {
  const auto &vertexA = vertices[facet.vertexA];
  float distance = facet.normal.x * (point.x - vertexA.x);
  distance += facet.normal.y * (point.y - vertexA.y);
  distance += facet.normal.z * (point.z - vertexA.z);
  return distance;
}
} // namespace

void Hull::update(const Coordinate &vertex_of_new_cone) {
  // find first visible facet
  std::size_t pos = 0;
  for (const auto &facet : vertices_and_faces.faces) {
    if (facet_point_distance(vertices_and_faces.vertices, *facet,
                             vertex_of_new_cone) > HULL_GEOMETRIC_TOLLERANCE) {
      this->update_(vertex_of_new_cone, pos);
      return;
    }
    ++pos;
  }
  throw Error{"Vertex passed to update hull is not outside it"};
}

void Hull::update(const Coordinate &vertex_of_new_cone,
                  const std::size_t starting_facet_for_expansion) {
  const auto &faces = vertices_and_faces.faces;
  if (faces.size() <= starting_facet_for_expansion) {
    throw Error{"Out of bounds facet index"};
  }
  if (facet_point_distance(vertices_and_faces.vertices,
                           *faces[starting_facet_for_expansion],
                           vertex_of_new_cone) <= HULL_GEOMETRIC_TOLLERANCE) {
    throw Error{"The specified starting facet is not valid"};
  }
  this->update_(vertex_of_new_cone, starting_facet_for_expansion);
}

namespace {
struct Edge {
  std::size_t vertex_first;
  std::size_t vertex_second;

  Facet *non_visible_neighbour_face;
  enum ConnectivityCase { AB, BC, CA };
  ConnectivityCase non_visible_neighbour_connectivity_case;
};

Edge::ConnectivityCase find_connectivity_case(const Facet &subject,
                                              const Facet *neighbour_to_find) {
  if (neighbour_to_find == subject.neighbourAB) {
    return Edge::ConnectivityCase::AB;
  }
  if (neighbour_to_find == subject.neighbourBC) {
    return Edge::ConnectivityCase::BC;
  }
  if (neighbour_to_find == subject.neighbourCA) {
    return Edge::ConnectivityCase::CA;
  }
  throw Error{"neighbour index not found"};
}
} // namespace

struct Hull::VisibleCone {
  std::vector<Edge> edges;
  std::set<Facet *> visible_faces;
};

Hull::VisibleCone Hull::computeVisibleCone(const Coordinate &vertex_of_new_cone,
                                           const std::size_t starting_facet) {
  auto &facets = vertices_and_faces.faces;
  const auto &vertices = vertices_and_faces.vertices;

  std::set<Facet *> visible_group = {};
  std::list<Facet *> facets_to_vist = {
      facets[starting_facet]
          .get()}; // this set contains facets we know that are visible, but
                   // whose neighbouring was not already visited
  std::vector<Edge> edges;
  while (!facets_to_vist.empty()) {
    auto *to_visit = facets_to_vist.front();
    facets_to_vist.pop_front();
    if (visible_group.find(to_visit) != visible_group.end()) {
      continue;
    }
    visible_group.emplace(to_visit);

    {
      auto *neighbourAB = to_visit->neighbourAB;
      if (facet_point_distance(vertices, *neighbourAB, vertex_of_new_cone) >
          HULL_GEOMETRIC_TOLLERANCE) {
        facets_to_vist.push_back(neighbourAB);
      } else {
        edges.push_back(Edge{to_visit->vertexA, to_visit->vertexB, neighbourAB,
                             find_connectivity_case(*neighbourAB, to_visit)});
      }
    }

    {
      auto *neighbourBC = to_visit->neighbourBC;
      if (facet_point_distance(vertices, *neighbourBC, vertex_of_new_cone) >
          HULL_GEOMETRIC_TOLLERANCE) {
        facets_to_vist.push_back(neighbourBC);
      } else {
        edges.push_back(Edge{to_visit->vertexB, to_visit->vertexC, neighbourBC,
                             find_connectivity_case(*neighbourBC, to_visit)});
      }
    }

    {
      auto *neighbourCA = to_visit->neighbourCA;
      if (facet_point_distance(vertices, *neighbourCA, vertex_of_new_cone) >
          HULL_GEOMETRIC_TOLLERANCE) {
        facets_to_vist.push_back(neighbourCA);
      } else {
        edges.push_back(Edge{to_visit->vertexC, to_visit->vertexA, neighbourCA,
                             find_connectivity_case(*neighbourCA, to_visit)});
      }
    }
  }
  return VisibleCone{std::move(edges), std::move(visible_group)};
}

namespace {
std::size_t find_edge_sharing_vertex(const std::vector<Edge> &edges,
                                     const std::size_t vertex,
                                     const std::size_t index_to_skip) {
  std::size_t result = 0;
  for (const auto &edge : edges) {
    if (((edge.vertex_first == vertex) || (edge.vertex_second == vertex)) &&
        (result != index_to_skip)) {
      return result;
    }
    ++result;
  }
  throw Error{"Vertex not found"};
}
} // namespace

void Hull::update_(const Coordinate &vertex_of_new_cone,
                   const std::size_t starting_facet_for_expansion) {
  auto visibility_cone =
      computeVisibleCone(vertex_of_new_cone, starting_facet_for_expansion);
  std::vector<Facet *> changed_facets;
  std::vector<Facet *> added_facets;
  RemovedFacets removed_facets;

  int facets_to_add = static_cast<int>(visibility_cone.edges.size()) -
                      static_cast<int>(visibility_cone.visible_faces.size());
  changed_facets = std::vector<Facet *>{visibility_cone.visible_faces.begin(),
                                        visibility_cone.visible_faces.end()};
  if (facets_to_add >= 0) {
    added_facets.reserve(facets_to_add);
    for (int k = 0; k < facets_to_add; ++k) {
      auto *added_facet =
          vertices_and_faces.faces.emplace_back(std::make_unique<Facet>())
              .get();
      added_facets.push_back(added_facet);
    }
  } else {
    const std::size_t remaining_facets = changed_facets.size() + facets_to_add;
    auto it_changed = changed_facets.begin();
    std::advance(it_changed, remaining_facets);

    removed_facets.reserve(-facets_to_add);
    std::for_each(it_changed, changed_facets.end(),
                  [&removed_facets,
                   &facets = vertices_and_faces.faces](hull::Facet *element) {
                    FacetPtr removed_facet = std::make_unique<Facet>(*element);
                    removed_facet->neighbourAB = nullptr;
                    removed_facet->neighbourBC = nullptr;
                    removed_facet->neighbourCA = nullptr;
                    removed_facets.emplace_back(std::move(removed_facet));
                    auto facets_it =
                        std::find_if(facets.begin(), facets.end(),
                                     [&element](const FacetPtr &facet_element) {
                                       return facet_element.get() == element;
                                     });
                    facets.erase(facets_it);
                  });

    changed_facets = std::vector<Facet *>{changed_facets.begin(), it_changed};
  }

  auto cone_facets = changed_facets;
  cone_facets.reserve(changed_facets.size() + added_facets.size());
  for (auto *facet : added_facets) {
    cone_facets.push_back(facet);
  }

  // build cone of new facets
  const std::size_t new_vertex_index = vertices_and_faces.vertices.size();
  vertices_and_faces.vertices.push_back(vertex_of_new_cone);
  std::size_t edge_index = 0;
  for (const auto &edge : visibility_cone.edges) {
    auto *facet_to_build = cone_facets[edge_index];
    facet_to_build->vertexA = edge.vertex_first;
    facet_to_build->vertexB = edge.vertex_second;
    facet_to_build->vertexC = new_vertex_index;

    facet_to_build->neighbourAB = edge.non_visible_neighbour_face;
    switch (edge.non_visible_neighbour_connectivity_case) {
    case Edge::ConnectivityCase::AB:
      facet_to_build->neighbourAB->neighbourAB = facet_to_build;
      break;
    case Edge::ConnectivityCase::BC:
      facet_to_build->neighbourAB->neighbourBC = facet_to_build;
      break;
    case Edge::ConnectivityCase::CA:
      facet_to_build->neighbourAB->neighbourCA = facet_to_build;
      break;
    }

    facet_to_build->neighbourCA = cone_facets[find_edge_sharing_vertex(
        visibility_cone.edges, facet_to_build->vertexA, edge_index)];

    facet_to_build->neighbourBC = cone_facets[find_edge_sharing_vertex(
        visibility_cone.edges, facet_to_build->vertexB, edge_index)];

    recomputeNormal(*facet_to_build);
    ++edge_index;
  }

  // update mid point for better computations
  const auto &vertices = vertices_and_faces.vertices;
  prod(Mid_point, static_cast<float>(vertices.size() - 1));
  Mid_point.x += vertices.back().x;
  Mid_point.y += vertices.back().y;
  Mid_point.z += vertices.back().z;
  prod(Mid_point, 1.f / static_cast<float>(vertices.size()));

  if (nullptr != this->observer) {
    observer->hullChanges(Observer::Notification{
        std::vector<const Facet *>{changed_facets.begin(),
                                   changed_facets.end()},
        std::vector<const Facet *>{added_facets.begin(), added_facets.end()},
        std::move(removed_facets), vertices_and_faces});
  }
}
} // namespace hull
