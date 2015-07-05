cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "-std=c++11 -D__TESTING__ -pthread")

project(GA_TSP_TEST)

file(GLOB_RECURSE TEST_SOURCES "cpp/*.cpp" "../src/*.cpp")
file(GLOB_RECURSE TEST_HEADERS "hpp/*.hpp" "../inc/*.hpp")

set (TEST_INCLUDE_DIRS "")
foreach (_headerFile ${TEST_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND TEST_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES TEST_INCLUDE_DIRS)

include_directories(${TEST_INCLUDE_DIRS})
add_executable (test testerMain.cpp ${TEST_SOURCES})

target_link_libraries(test pthread boost_system boost_serialization boost_random)
