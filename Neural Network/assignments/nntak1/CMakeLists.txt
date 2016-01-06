cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "-std=c++11 -pthread -DDEBUG")

project(NNTAK1)

file(GLOB_RECURSE ISOURCES "src/*.cpp")
file(GLOB_RECURSE IHEADER "inc/*.hxx")

set (IINCLUDE_DIRS "")
foreach (_headerFile ${IHEADER})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND IINCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES IINCLUDE_DIRS)

include_directories(${IINCLUDE_DIRS} "." "vendors/eigen/")

add_executable (nn main.cpp ${ISOURCES})

target_link_libraries (nn pthread boost_system boost_program_options boost_filesystem matio)

# set warnings level
IF(MSVC) # MSVC compiler (Win32 only)
    # Display more warnings
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
ELSEIF(UNIX OR CMAKE_COMPILER_IS_GNUCXX) # Clang OR Gcc (Linux, Mac OS or Win32 with MingW)
    # displays all warnings
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
ENDIF()