set(PROJECT_SHORTNAME Hull)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(${PROJECT_SHORTNAME} ${SRC_FILES})

target_include_directories(${PROJECT_SHORTNAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/header/)
