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

void Hull::setObserver(Observer &obs) { this->observer = &obs; }

void Hull::recomputeNormal(Facet &subject) const {
  Coordinate delta1, delta2;
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
        vertices});
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
  for (const auto &facet : facets) {
    if (facet_point_distance(vertices, *facet, vertex_of_new_cone) >
        HULL_GEOMETRIC_TOLLERANCE) {
      this->update_(vertex_of_new_cone, pos);
      return;
    }
    ++pos;
  }
  throw Error{"Vertex passed to update hull is not outside it"};
}

void Hull::update(const Coordinate &vertex_of_new_cone,
                  const std::size_t starting_facet_for_expansion) {
  if (facets.size() <= starting_facet_for_expansion) {
    throw Error{"Out of bounds facet index"};
  }
  if (facet_point_distance(vertices, *facets[starting_facet_for_expansion],
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

Edge::ConnectivityCase find_connectivity_case(Facet &subject,
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
  throw std::runtime_error{"neighbour index not found"};
}
} // namespace

struct Hull::VisibleCone {
  std::vector<Edge> edges;
  std::vector<std::size_t> visible_faces;
};

Hull::VisibleCone Hull::computeVisibleCone(const Coordinate &vertex_of_new_cone,
                                           const std::size_t starting_facet) {
  std::set<std::size_t> visible_group = {};
  std::list<std::size_t> open_set = {
      starting_facet}; // this set contains facets we know that are visible, but
                       // whose neighbouring was not already visited
  std::vector<Edge> edges;
  while (!open_set.empty()) {
    std::size_t to_visit = open_set.front();
    open_set.pop_front();
    if (visible_group.find(to_visit) != visible_group.end()) {
      continue;
    }
    visible_group.emplace(to_visit);

    const auto *neighbour_facet_AB = facets[to_visit]->neighbourAB;
    if (facet_point_distance(vertices, *neighbour_facet_AB,
                             vertex_of_new_cone) > HULL_GEOMETRIC_TOLLERANCE) {
      open_set.push_back(neighbourAB_index);
    } else {
      edges.push_back(Edge{
          facets[to_visit].vertexA, facets[to_visit].vertexB, neighbourAB_index,
          find_connectivity_case(facets[neighbourAB_index], to_visit)});
    }

    const auto &neighbourBC_index = facets[to_visit].neighbourBC;
    if (facet_point_distance(vertices, facets[neighbourBC_index],
                             vertex_of_new_cone) > HULL_GEOMETRIC_TOLLERANCE) {
      open_set.push_back(neighbourBC_index);
    } else {
      edges.push_back(Edge{
          facets[to_visit].vertexB, facets[to_visit].vertexC, neighbourBC_index,
          find_connectivity_case(facets[neighbourBC_index], to_visit)});
    }

    const auto &neighbourCA_index = facets[to_visit].neighbourCA;
    if (facet_point_distance(vertices, facets[neighbourCA_index],
                             vertex_of_new_cone) > HULL_GEOMETRIC_TOLLERANCE) {
      open_set.push_back(neighbourCA_index);
    } else {
      edges.push_back(Edge{
          facets[to_visit].vertexC, facets[to_visit].vertexA, neighbourCA_index,
          find_connectivity_case(facets[neighbourCA_index], to_visit)});
    }
  }
  return VisibleCone{
      std::move(edges),
      std::vector<std::size_t>{visible_group.begin(), visible_group.end()}};
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

void remove_facets(std::vector<Facet> &facets,
                   std::list<std::size_t> indices_to_remove) {
  while (!indices_to_remove.empty()) {
    std::size_t index = indices_to_remove.front();
    indices_to_remove.pop_front();
    auto facets_it = facets.begin();
    std::advance(facets_it, index);
    facets.erase(facets_it);
    for (auto &other_index : indices_to_remove) {
      if (other_index > index) {
        --other_index;
      }
    }
  }
}
} // namespace

void Hull::update_(const Coordinate &vertex_of_new_cone,
                   const std::size_t starting_facet_for_expansion) {
  auto visibility_cone =
      computeVisibleCone(vertex_of_new_cone, starting_facet_for_expansion);
  std::vector<std::size_t> changed_facets;
  std::vector<std::size_t> added_facets;
  std::vector<Observer::Notification::FacetAndOldIndex> removed_facets;

  int facets_to_add = static_cast<int>(visibility_cone.edges.size()) -
                      static_cast<int>(visibility_cone.visible_faces.size());
  changed_facets = visibility_cone.visible_faces;
  if (facets_to_add >= 0) {
    added_facets.reserve(facets_to_add);
    for (int k = 0; k < facets_to_add; ++k) {
      added_facets.push_back(facets.size());
      facets.emplace_back();
    }
  } else {
    std::size_t remaining_facets = changed_facets.size() + facets_to_add;
    auto it_changed = changed_facets.begin();
    std::advance(it_changed, remaining_facets);

    std::list<std::size_t> indices_to_remove(it_changed, changed_facets.end());
    removed_facets.reserve(indices_to_remove.size());
    for (const auto &index : indices_to_remove) {
      removed_facets.push_back(
          Observer::Notification::FacetAndOldIndex{facets[index], index});
    }
    remove_facets(facets, std::move(indices_to_remove));

    changed_facets =
        std::vector<std::size_t>{changed_facets.begin(), it_changed};
  }

  std::vector<std::size_t> cone_facets = changed_facets;
  cone_facets.reserve(changed_facets.size() + added_facets.size());
  for (const auto index : added_facets) {
    cone_facets.push_back(index);
  }

  // build cone of new facets
  const std::size_t new_vertex_index = vertices.size();
  vertices.push_back(vertex_of_new_cone);
  std::size_t edge_index = 0;
  for (const auto &edge : visibility_cone.edges) {
    std::size_t facet_to_build_index = cone_facets[edge_index];
    auto &facet_to_build = facets[facet_to_build_index];
    facet_to_build.vertexA = edge.vertex_first;
    facet_to_build.vertexB = edge.vertex_second;
    facet_to_build.vertexC = new_vertex_index;

    facet_to_build.neighbourAB = edge.non_visible_neighbour_face_index;
    switch (edge.non_visible_neighbour_connectivity_case) {
    case Edge::ConnectivityCase::AB:
      facets[edge.non_visible_neighbour_face_index].neighbourAB =
          facet_to_build_index;
      break;
    case Edge::ConnectivityCase::BC:
      facets[edge.non_visible_neighbour_face_index].neighbourBC =
          facet_to_build_index;
      break;
    case Edge::ConnectivityCase::CA:
      facets[edge.non_visible_neighbour_face_index].neighbourCA =
          facet_to_build_index;
      break;
    }

    facet_to_build.neighbourCA = cone_facets[find_edge_sharing_vertex(
        visibility_cone.edges, facet_to_build.vertexA, edge_index)];

    facet_to_build.neighbourBC = cone_facets[find_edge_sharing_vertex(
        visibility_cone.edges, facet_to_build.vertexB, edge_index)];

    recomputeNormal(facet_to_build);
    ++edge_index;
  }

  // update mid point for better computations
  prod(Mid_point, static_cast<float>(vertices.size() - 1));
  Mid_point.x += vertices.back().x;
  Mid_point.y += vertices.back().y;
  Mid_point.z += vertices.back().z;
  prod(Mid_point, 1.f / static_cast<float>(vertices.size()));

  if (nullptr != this->observer) {
    observer->hullChanges(Observer::Notification{
        std::move(changed_facets), std::move(added_facets),
        std::move(removed_facets), vertices, facets});
  }
}
} // namespace hull
