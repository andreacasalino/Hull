cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(HullUtils)
add_subdirectory(src)

option(BUILD_HULL_TESTS "" OFF)
if(BUILD_HULL_TESTS)
	project(Hull-Tests)
    add_subdirectory(tests)
endif()

