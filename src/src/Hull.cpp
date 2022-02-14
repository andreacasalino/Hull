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
    std::vector<const Facet *> initial_facets;
    initial_facets.reserve(4);
    for (const auto &facet : facets) {
      initial_facets.push_back(&facet);
    }
    observer->hullChanges(Observer::Notification{initial_facets, {}, {}});
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

std::vector<bool>
Hull::computeVisibilityFlags(const Coordinate& vertex_of_new_cone,
    const std::size_t starting_facet) const {
    std::set<std::size_t> visible_group = { starting_facet };
    std::list<std::size_t> open_set = { facets[starting_facet].neighbourAB,
                                       facets[starting_facet].neighbourBC,
                                       facets[starting_facet].neighbourCA };
    while (!open_set.empty()) {
        std::size_t to_visit = open_set.front();
        open_set.pop_front();
        if ((visible_group.find(to_visit) == visible_group.end()) &&
            (facet_point_distance(vertices, facets[to_visit], vertex_of_new_cone) >
                HULL_GEOMETRIC_TOLLERANCE)) {
            visible_group.emplace(to_visit);
            open_set.push_back(facets[to_visit].neighbourAB);
            open_set.push_back(facets[to_visit].neighbourBC);
            open_set.push_back(facets[to_visit].neighbourCA);
        }
    }
    std::vector<bool> flags(facets.size(), false);
    for (const auto index : visible_group) {
        flags[index] = true;
    }
    return flags;
}

namespace {
    void Replace(Facet& involved_facet, const std::size_t old_neigh_index, const std::size_t new_neigh_index) {
        if (involved_facet.neighbourAB == old_neigh_index) {
            involved_facet.neighbourAB = new_neigh_index;
        }
        else if (involved_facet.neighbourBC == old_neigh_index) {
            involved_facet.neighbourBC = new_neigh_index;
        }
        else {
            involved_facet.neighbourCA = new_neigh_index;
        }
    };
}

void Hull::update_(const Coordinate& vertex_of_new_cone, const std::size_t starting_facet_for_expansion) {
    auto visibility_flags = computeVisibilityFlags(vertex_of_new_cone, starting_facet_for_expansion);

}

void Hull::_UpdateHull(const Coordinate &vertex_of_new_cone,
                       Facet &starting_facet_for_expansion) {
  starting_facet_for_expansion.bVisible = true;
  const Coordinate *Point = &this->vertices.emplace_back(vertex_of_new_cone);
  // find the group of visible faces
  std::list<Facet *> Visible_group;
  Visible_group.push_back(&starting_facet_for_expansion);
  auto itN = Visible_group.begin();
  std::size_t kNeigh = 0, k;
  float distance;
  while (kNeigh < Visible_group.size()) {
    for (k = 0; k < 3; ++k) {
      if (!(*itN)->Neighbour[k]->bVisible) { // this neighbour facet is not
                                             // already present in Visible_Group
        distance =
            (*itN)->Neighbour[k]->N.x * (Point->x - (*itN)->Neighbour[k]->A->x);
        distance +=
            (*itN)->Neighbour[k]->N.y * (Point->y - (*itN)->Neighbour[k]->A->y);
        distance +=
            (*itN)->Neighbour[k]->N.z * (Point->z - (*itN)->Neighbour[k]->A->z);
        if (distance > HULL_GEOMETRIC_TOLLERANCE) {
          (*itN)->Neighbour[k]->bVisible = true;
          Visible_group.push_back((*itN)->Neighbour[k]);
        }
      }
    }
    ++kNeigh;
    ++itN;
  }
  // Update Hull. Build the cone of new facets
  const Coordinate *A, *B, *C;
  Facet *AB, *BC, *CA;
  itN = Visible_group.begin();
  std::size_t kboard, initialVisibleSize = Visible_group.size();
  std::list<const Facet *> removed, changed, created;
  for (k = 0; k < initialVisibleSize; ++k) {
    kboard = 0;
    A = (*itN)->A;
    B = (*itN)->B;
    C = (*itN)->C;
    AB = (*itN)->Neighbour[0];
    BC = (*itN)->Neighbour[1];
    CA = (*itN)->Neighbour[2];
    // AB
    if (!AB->bVisible) { // edge AB is part of the board of the cone
      (*itN)->C = Point;
      this->RecomputeNormal(**itN);
      ++kboard;
    }
    // BC
    if (!BC->bVisible) { // edge BC is part of the board of the cone
      if (kboard == 0) {
        (*itN)->A = B;
        (*itN)->B = C;
        (*itN)->C = Point;
        this->RecomputeNormal(**itN);
        (*itN)->Neighbour[0] = BC;
      } else {
        this->AppendFacet(*B, *C, *Point);
        this->Facets.back().bVisible = true;
        Visible_group.push_back(&this->Facets.back());
        created.push_back(&this->Facets.back());
        Visible_group.back()->Neighbour[0] = BC;
        Replace(*BC, **itN, *Visible_group.back());
      }
      ++kboard;
    }
    // CA
    if (!CA->bVisible) { // edge CA is part of the board of the cone
      if (kboard == 0) {
        (*itN)->A = A;
        (*itN)->B = C;
        (*itN)->C = Point;
        this->RecomputeNormal(**itN);
        (*itN)->Neighbour[0] = CA;
      } else {
        this->AppendFacet(*A, *C, *Point);
        this->Facets.back().bVisible = true;
        Visible_group.push_back(&this->Facets.back());
        created.push_back(&this->Facets.back());
        Visible_group.back()->Neighbour[0] = CA;
        Replace(*CA, **itN, *Visible_group.back());
      }
      ++kboard;
    }
    if (kboard == 0) {
      (*itN)->Neighbour[0] = nullptr;
      removed.push_back(*itN);
      itN = Visible_group.erase(itN);
    } else {
      changed.push_back(*itN);
      ++itN;
    }
  }
  if (nullptr != this->observer) {
    // notify to observer
    this->observer->RemovedFacets(removed);
  }
  // remove facet for which Neighbour[0]=nullptr
  auto itF = this->Facets.begin();
  while (itF != this->Facets.end()) {
    if (nullptr == itF->Neighbour[0])
      itF = this->Facets.erase(itF);
    else
      ++itF;
  }
  // delete old neighbour info for the visible group
  for (itN = Visible_group.begin(); itN != Visible_group.end(); ++itN) {
    (*itN)->Neighbour[1] = nullptr;
    (*itN)->Neighbour[2] = nullptr;
  }
  // update neighbour info for the visible group
  auto itN2 = itN;
  for (itN = Visible_group.begin(); itN != Visible_group.end(); ++itN) {
    if (nullptr == (*itN)->Neighbour[1]) { // find neighbour of edge BC
      for (itN2 = Visible_group.begin(); itN2 != Visible_group.end(); ++itN2) {
        if (*itN2 != *itN) {
          if ((*itN2)->A == (*itN)->B) {
            (*itN)->Neighbour[1] = *itN2;
            (*itN2)->Neighbour[2] = *itN;
            break;
          }

          if ((*itN2)->B == (*itN)->B) {
            (*itN)->Neighbour[1] = *itN2;
            (*itN2)->Neighbour[1] = *itN;
            break;
          }
        }
      }
    }

    if (nullptr == (*itN)->Neighbour[2]) { // find neighbour of edge CA
      for (itN2 = Visible_group.begin(); itN2 != Visible_group.end(); ++itN2) {
        if (*itN2 != *itN) {
          if ((*itN2)->A == (*itN)->A) {
            (*itN)->Neighbour[2] = *itN2;
            (*itN2)->Neighbour[2] = *itN;
            break;
          }

          if ((*itN2)->B == (*itN)->A) {
            (*itN)->Neighbour[2] = *itN2;
            (*itN2)->Neighbour[1] = *itN;
            break;
          }
        }
      }
    }
  }
  // clean up
  for (itN = Visible_group.begin(); itN != Visible_group.end(); ++itN)
    (*itN)->bVisible = false;
  if (nullptr != this->observer) {
    // notify to observer
    this->observer->AddedChangedFacets(created, changed);
  }
}
} // namespace hull
