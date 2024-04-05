/**
 * Author:    Andrea Casalino
 * Created:   03.12.2019
 *
 * report any bug to andrecasa91@gmail.com.
 **/

#pragma once

#include <Hull/Coordinate.h>

#include <deque>
#include <list>
#include <unordered_set>
#include <vector>

namespace hull {
struct Facet {
  std::size_t vertexA;
  std::size_t vertexB;
  std::size_t vertexC;

  Facet *neighbourAB;
  Facet *neighbourBC;
  Facet *neighbourCA;

  Coordinate normal; // outer normal
};

struct HullContext {
  std::deque<Coordinate> vertices;
  std::unordered_set<Facet *> faces;
  Coordinate internalPoint;
};

class Observer {
public:
  virtual ~Observer() = default;

  struct Notification {
    Notification(const HullContext &context) : context{context} {}

    std::vector<const Facet *> changed;
    std::vector<const Facet *> added;
    std::vector<const Facet *> removed;

    const HullContext &context;
  };

  virtual void hullChanges(const Notification &notification) = 0;
};

class FacetAllocator {
public:
  FacetAllocator();

  Facet *getFacet();

  void markAsRecyclable(Facet *facet) { outstanding.push_back(facet); }

private:
  struct Buffer {
    Buffer(std::size_t capacity);
    ~Buffer();

    std::size_t size = 0;
    const std::size_t capacity;
    char *buffer;
  };
  std::list<Buffer> buffers;
  std::deque<Facet *> outstanding;
};

class Hull {
public:
  Hull(const Coordinate &A, const Coordinate &B, const Coordinate &C,
       const Coordinate &D);

  Hull(const Coordinate &A, const Coordinate &B, const Coordinate &C,
       const Coordinate &D, Observer &obs);

  void setObserver(Observer &obs);

  void update(const Coordinate &vertex_of_new_cone);

  void update(const Coordinate &vertex_of_new_cone,
              Facet *starting_facet_for_expansion);

  const auto &getContext() const { return context; }

private:
  void initThetraedron(const Coordinate &A, const Coordinate &B,
                       const Coordinate &C, const Coordinate &D);

  void update_(const Coordinate &vertex_of_new_cone,
               Facet *starting_facet_for_expansion);

  Facet *makeFacet(std::size_t vertexA, std::size_t vertexB,
                   std::size_t vertexC);

  HullContext context;
  FacetAllocator allocator;

  Observer *observer = nullptr;
};
} // namespace hull
