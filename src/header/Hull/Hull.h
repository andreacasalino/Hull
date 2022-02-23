/**
 * Author:    Andrea Casalino
 * Created:   03.12.2019
 *
 * report any bug to andrecasa91@gmail.com.
 **/

#pragma once

#include <Hull/Coordinate.h>
#include <cmath>
#include <memory>
#include <vector>

namespace hull {
struct Facet {
  std::size_t vertexA;
  std::size_t vertexB;
  std::size_t vertexC;

  Facet *neighbourAB = nullptr;
  Facet *neighbourBC = nullptr;
  Facet *neighbourCA = nullptr;

  Coordinate normal; // outer normal
};
using FacetPtr = std::unique_ptr<Facet>;

using Facets = std::vector<FacetPtr>;
using RemovedFacets = std::vector<std::unique_ptr<const Facet>>;

class Observer {
public:
  virtual ~Observer() = default;

  struct Notification {
    std::vector<const Facet *> changed;

    std::vector<const Facet *> added;
    RemovedFacets removed;

    // context
    const std::vector<Coordinate> &vertices;
  };

  virtual void hullChanges(const Notification &notification) = 0;
};

class Hull {
public:
  Hull(const Coordinate &A, const Coordinate &B, const Coordinate &C,
       const Coordinate &D);

  void setObserver(Observer &obs);

  void update(const Coordinate &vertex_of_new_cone);

  void update(const Coordinate &vertex_of_new_cone,
              const std::size_t starting_facet_for_expansion);

  const std::vector<Coordinate> &getVertices() const { return this->vertices; };
  const Facets &getFacets() const { return this->facets; };

private:
  void initThetraedron(const Coordinate &A, const Coordinate &B,
                       const Coordinate &C, const Coordinate &D);

  FacetPtr Hull::makeFacet(const std::size_t vertexA, const std::size_t vertexB,
                           const std::size_t vertexC) const;

  void update_(const Coordinate &vertex_of_new_cone,
               const std::size_t starting_facet_for_expansion);

  void recomputeNormal(Facet &subject) const;

  FacetPtr makeFacet(const std::size_t vertexA, const std::size_t vertexB,
                     const std::size_t vertexC) const;

  struct VisibleCone;
  VisibleCone computeVisibleCone(const Coordinate &vertex_of_new_cone,
                                 const std::size_t starting_facet);

  std::vector<Coordinate> vertices;
  Facets facets;
  Coordinate Mid_point;
  Observer *observer = nullptr;
};
} // namespace hull
