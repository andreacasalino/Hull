/**
 * Author:    Andrea Casalino
 * Created:   03.12.2019
 *
 * report any bug to andrecasa91@gmail.com.
 **/

#include <Hull/Definitions.h>
#include <Hull/Error.h>
#include <Hull/Hull.h>
#include <limits>
#include <list>
#include <set>
#include <algorithm>

namespace hull {
Hull::Hull(const Coordinate &A, const Coordinate &B, const Coordinate &C,
           const Coordinate &D) {
  this->initThetraedron(A, B, C, D);
};

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

namespace {
constexpr std::size_t INVALID_FACET_INDEX =
    std::numeric_limits<std::size_t>::max();
}

Facet Hull::makeFacet(const std::size_t vertexA, const std::size_t vertexB,
                      const std::size_t vertexC) const {
  Facet result;
  result.vertexA = vertexA;
  result.vertexB = vertexB;
  result.vertexC = vertexC;
  recomputeNormal(result);
  result.neighbourAB = INVALID_FACET_INDEX;
  result.neighbourBC = INVALID_FACET_INDEX;
  result.neighbourCA = INVALID_FACET_INDEX;
  return result;
}

void Hull::initThetraedron(const Coordinate &A, const Coordinate &B,
                           const Coordinate &C, const Coordinate &D) {
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

  // computation of the midpoint of the thetraedron
  this->Mid_point.x = 0.25f * (A.x + B.x + C.x + D.x);
  this->Mid_point.y = 0.25f * (A.y + B.y + C.y + D.y);
  this->Mid_point.z = 0.25f * (A.z + B.z + C.z + D.z);

  // build the tethraedron
  // ABC->0; ABD->1; ACD->2; BCD->3
  facets.reserve(4);
  facets.push_back(makeFacet(0, 1, 2));
  facets.push_back(makeFacet(0, 1, 3));
  facets.push_back(makeFacet(0, 2, 3));
  facets.push_back(makeFacet(1, 2, 3));

  // setup initial connectivity
  // ABC
  facets[0].neighbourAB = 1;
  facets[0].neighbourBC = 3;
  facets[0].neighbourCA = 2;
  // ABD
  facets[1].neighbourAB = 0;
  facets[1].neighbourBC = 3;
  facets[1].neighbourCA = 2;
  // ACD
  facets[2].neighbourAB = 0;
  facets[2].neighbourBC = 3;
  facets[2].neighbourCA = 1;
  // BCD
  facets[3].neighbourAB = 0;
  facets[3].neighbourBC = 2;
  facets[3].neighbourCA = 1;

  if (nullptr != this->observer) {
    observer->hullChanges(Observer::Notification{ {0,1,2,3}, {}, {} });
  }
}

namespace {
float facet_point_distance(const std::vector<Coordinate> &vertices,
                           const Facet &facet, const Coordinate &point) {
  const auto &vertexACoord = vertices[facet.vertexA];
  float distance = facet.normal.x * (point.x - vertexACoord.x);
  distance += facet.normal.y * (point.y - vertexACoord.y);
  distance += facet.normal.z * (point.z - vertexACoord.z);
  return distance;
}
} // namespace

void Hull::update(const Coordinate& vertex_of_new_cone) {
    // find first visible facet
    std::size_t pos = 0;
    for (const auto& facet : facets) {
        if (facet_point_distance(vertices, facet, vertex_of_new_cone) > HULL_GEOMETRIC_TOLLERANCE) {
            this->update_(vertex_of_new_cone, pos);
            return;
        }
        ++pos;
    }
    throw Error{ "Vertex passed to update hull is not outside it" };
}

void Hull::update(const Coordinate& vertex_of_new_cone, const std::size_t starting_facet_for_expansion) {
    if (facets.size() <= starting_facet_for_expansion) {
        throw Error{ "Out of bounds facet index" };
    }
    if (facet_point_distance(vertices, facets[starting_facet_for_expansion], vertex_of_new_cone) <= HULL_GEOMETRIC_TOLLERANCE) {
        throw Error{ "The specified starting facet is not valid" };
    }
    this->update_(vertex_of_new_cone, starting_facet_for_expansion);
}

Hull::VisibleCone Hull::computeVisibleCone(const Coordinate& vertex_of_new_cone, const std::size_t starting_facet) const {
    std::set<std::size_t> visible_group = {};
    std::list<std::size_t> open_set = { starting_facet }; // this set contain facet we know are visible, but whose neighbouring was not already computed
    std::vector<Edge> edges;
    while (!open_set.empty()) {
        std::size_t to_visit = open_set.front();
        open_set.pop_front();
        if (visible_group.find(to_visit) != visible_group.end()) {
            continue;
        }
        visible_group.emplace(to_visit);

        const auto& to_visit_neighbourAB = facets[to_visit].neighbourAB;
        if (facet_point_distance(vertices, facets[to_visit_neighbourAB], vertex_of_new_cone) > HULL_GEOMETRIC_TOLLERANCE) {
            open_set.push_back(to_visit_neighbourAB);
        }
        else {
            edges.push_back(Edge{facets[to_visit].vertexA, facets[to_visit].vertexB, to_visit_neighbourAB });
        }

        const auto& to_visit_neighbourBC = facets[to_visit].neighbourBC;
        if (facet_point_distance(vertices, facets[to_visit_neighbourBC], vertex_of_new_cone) > HULL_GEOMETRIC_TOLLERANCE) {
            open_set.push_back(to_visit_neighbourBC);
        }
        else {
            edges.push_back(Edge{ facets[to_visit].vertexB, facets[to_visit].vertexC, to_visit_neighbourBC });
        }

        const auto& to_visit_neighbourCA = facets[to_visit].neighbourCA;
        if (facet_point_distance(vertices, facets[to_visit_neighbourCA], vertex_of_new_cone) > HULL_GEOMETRIC_TOLLERANCE) {
            open_set.push_back(to_visit_neighbourCA);
        }
        else {
            edges.push_back(Edge{ facets[to_visit].vertexC, facets[to_visit].vertexA, to_visit_neighbourCA });
        }
    }
    return VisibleCone{std::move(edges), std::vector<std::size_t>{visible_group.begin(), visible_group.end()}};
}

namespace {
    std::size_t find_edge_sharing_vertex(const std::vector<Hull::Edge>& edges, const std::size_t vertex, const Hull::Edge& involved_edge) {
        std::size_t result = 0;
        for (const auto& edge: edges) {
            if ((edge.vertex_first == vertex) || (edge.vertex_second == vertex)
                && (&edge != &involved_edge)) {
                return result;
            }
            ++result;
        }
        throw Error{"Vertex not found"};
    }
}

void Hull::update_(const Coordinate& vertex_of_new_cone, const std::size_t starting_facet_for_expansion) {
    auto visibility_cone = computeVisibleCone(vertex_of_new_cone, starting_facet_for_expansion);
    std::vector<std::size_t> changed_facets;
    std::vector<std::size_t> added_facets;
    std::vector<Facet> removed_facets;

    int facets_to_add = static_cast<int>(visibility_cone.edges.size()) - static_cast<int>(visibility_cone.visible_faces.size());
    changed_facets = visibility_cone.visible_faces;
    if (facets_to_add >= 0) {
        added_facets.reserve(facets_to_add);
        for (std::size_t k = 0; k < facets_to_add; ++k) {
            added_facets.push_back(facets.size());
            facets.emplace_back();
        }
    }
    else {
        std::size_t remaining_facets = changed_facets.size() + facets_to_add;
        auto it_changed = changed_facets.begin();
        std::advance(it_changed, remaining_facets);
        removed_facets.reserve(-facets_to_add);
        std::for_each(it_changed, changed_facets.end(), [&facets = this->facets, &removed_facets](const std::size_t facet_pos) {
            auto it_facets = facets.begin();
            std::advance(it_facets, facet_pos);
            removed_facets.push_back(*it_facets);
            facets.erase(it_facets);
        });
        changed_facets = std::vector<std::size_t>{ changed_facets.begin(), it_changed};
    }

    std::vector<std::size_t> all_remaining_facets = changed_facets;
    all_remaining_facets.reserve(changed_facets.size() + added_facets.size());
    for (const auto index : added_facets) {
        all_remaining_facets.push_back(index);
    }

    // build cone of new facets
    const std::size_t new_vertex_index = vertices.size();
    vertices.push_back(vertex_of_new_cone);
    std::size_t emplacing_index = 0;
    for (const auto& edge: visibility_cone.edges) {
        auto& facet_to_build = facets[all_remaining_facets[emplacing_index]];
        facet_to_build.vertexA = edge.vertex_first;
        facet_to_build.vertexB = edge.vertex_second;
        facet_to_build.vertexC = new_vertex_index;

        facet_to_build.neighbourAB = edge.neighbour_face;
        facet_to_build.neighbourCA = all_remaining_facets[find_edge_sharing_vertex(visibility_cone.edges, facet_to_build.vertexA, edge)];
        facet_to_build.neighbourBC = all_remaining_facets[find_edge_sharing_vertex(visibility_cone.edges, facet_to_build.vertexB, edge)];

        recomputeNormal(facet_to_build);
        ++emplacing_index;
    }

    if (nullptr != this->observer) {
        observer->hullChanges(Observer::Notification{ std::move(added_facets), std::move(changed_facets), std::move(removed_facets) });
    }
}
} // namespace hull
