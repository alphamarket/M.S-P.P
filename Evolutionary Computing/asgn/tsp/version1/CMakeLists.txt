cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "-std=c++11 -D__DEBUG__ -pthread")

project(GA_TSP)

file(GLOB_RECURSE TSP_SOURCES "src/*.cpp")
file(GLOB_RECURSE TSP_HEADERS "inc/*.hpp")

set (TSP_INCLUDE_DIRS "")
foreach (_headerFile ${TSP_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND TSP_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES TSP_INCLUDE_DIRS)

include_directories(${TSP_INCLUDE_DIRS})
add_executable (main main.cpp ${TSP_SOURCES})

target_link_libraries (main
    boost_program_options pthread boost_regex boost_system boost_filesystem boost_serialization)
    
include_directories(vendors/jsoncons/src)
