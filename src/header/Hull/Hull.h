/**
 * Author:    Andrea Casalino
 * Created:   03.12.2019
 *
 * report any bug to andrecasa91@gmail.com.
 **/

#pragma once

#include <Hull/Coordinate.h>
#include <cmath>
#include <vector>

namespace hull {
struct Facet {
  std::size_t vertexA;
  std::size_t vertexB;
  std::size_t vertexC;

  std::size_t neighbourAB;
  std::size_t neighbourBC;
  std::size_t neighbourCA;

  Coordinate normal; // outer normal
};

class Observer {
public:
  virtual ~Observer() = default;

  struct Notification {
    std::vector<std::size_t> added;
    std::vector<std::size_t> changed;
    std::vector<Facet> removed;
  };

  virtual void hullChanges(const Notification &notification) = 0;
};

class Hull {
public:
  Hull(const Coordinate &A, const Coordinate &B, const Coordinate &C,
       const Coordinate &D);

  void setObserver(Observer &obs);

  void update(const Coordinate& vertex_of_new_cone);

  void update(const Coordinate &vertex_of_new_cone,
              const std::size_t starting_facet_for_expansion);

  const std::vector<Coordinate> &getVertices() const { return this->vertices; };
  const std::vector<Facet> &getFacets() const { return this->facets; };

  struct Edge {
      std::size_t vertex_first;
      std::size_t vertex_second;
      std::size_t neighbour_face;
  };
  struct VisibleCone {
      std::vector<Edge> edges;
      std::vector<std::size_t> visible_faces;
  };

private:
  void initThetraedron(const Coordinate &A, const Coordinate &B,
                       const Coordinate &C, const Coordinate &D);

  void update_(const Coordinate& vertex_of_new_cone, const std::size_t starting_facet_for_expansion);

  void recomputeNormal(Facet &subject) const;

  Facet makeFacet(const std::size_t vertexA, const std::size_t vertexB,
                  const std::size_t vertexC) const;

  VisibleCone computeVisibleCone(const Coordinate &vertex_of_new_cone, const std::size_t starting_facet) const;

  std::vector<Coordinate> vertices;
  std::vector<Facet> facets;
  Coordinate Mid_point;
  Observer *observer = nullptr;
};
} // namespace hull
