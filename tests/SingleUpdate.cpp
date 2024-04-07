#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "Utils.h"
#include <Hull/Coordinate.h>

TEST_CASE("Hull single udpate") {
  std::vector<hull::Coordinate> initial_vertices = {
      {1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, -1}};

  SECTION("Single visible facet") {
    hull::Coordinate vertex_for_update = {0.5, 0, 1};

    hull::TrivialObserver obs;
    hull::Hull hull(initial_vertices[0], initial_vertices[1],
                    initial_vertices[2], initial_vertices[3]);
    hull.setObserver(obs);

    hull.update(vertex_for_update);

    // check the update
    CHECK(obs.last_notification->changed.size() == 1);
    CHECK(obs.last_notification->added.size() == 2);
    CHECK(obs.last_notification->removed.empty());
    CHECK(hull.getContext().vertices.size() == 5);
    CHECK(hull.getContext().faces.size() == 6);
  }

  SECTION("Two visible facets") {
    hull::Coordinate vertex_for_update = {-0.5, 0, 1};

    hull::TrivialObserver obs;
    hull::Hull hull(initial_vertices[0], initial_vertices[1],
                    initial_vertices[2], initial_vertices[3]);
    hull.setObserver(obs);

    hull.update(vertex_for_update);

    // check the update
    CHECK(obs.last_notification->changed.size() == 2);
    CHECK(obs.last_notification->added.size() == 2);
    CHECK(obs.last_notification->removed.empty());
    CHECK(hull.getContext().vertices.size() == 5);
    CHECK(hull.getContext().faces.size() == 6);
  }
}
