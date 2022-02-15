#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <Hull/Coordinate.h>
#include "HullToObj.h"

#include <sstream>
#include <algorithm>

hull::Hull make_hull(const std::vector<hull::Coordinate>& vertices) {
	hull::Hull hull(vertices[0], vertices[1], vertices[12], vertices[3]);
	auto it_vertices = vertices.begin();
	std::advance(it_vertices, 4);
	std::for_each(it_vertices, vertices.end(), [&hull](const hull::Coordinate& vertex) {
		hull.update(vertex);
	});
	return hull;
}

TEST_CASE("Simple thetraedron") {

}

TEST_CASE("Invalid thetraedrons") {

}

std::vector<hull::Coordinate> make_sphere_cloud(const std::size_t samples);

TEST_CASE("Sphere clouds") {
	std::size_t angles_samples = GENERATE(5, 10, 20, 50, 100, 500);

	std::vector<hull::Coordinate> vertices = make_sphere_cloud(samples);

	auto hull = make_hull(vertices);

	CHECK(hull.getVertices().size() == samples);
	// check all the normals are pointing outside the origin
	for (const auto& facet : hull.getFacets()) {
		CHECK(hull::dot(hull.getVertices()[facet.vertexA], facet.normal) > 0);
	}

	std::stringstream log_name;
	log_name << "SphereCloud-" << std::to_string(samples) << ".obj";
	hull::toObj(hull, log_name.str());
}
