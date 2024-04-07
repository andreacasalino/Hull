include(FetchContent)
FetchContent_Declare(
	catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        v3.4.0
)
FetchContent_MakeAvailable(catch2)
	
set(TEST_NAME Tests)

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

add_executable(${TEST_NAME} ${SOURCES})

target_link_libraries(${TEST_NAME} PUBLIC
    Catch2::Catch2
    Catch2::Catch2WithMain
    Hull
)

target_compile_definitions(${TEST_NAME} PUBLIC
    -D LOG_FOLDER="${CMAKE_CURRENT_SOURCE_DIR}/LogFolder"
)

install(TARGETS ${TEST_NAME})
