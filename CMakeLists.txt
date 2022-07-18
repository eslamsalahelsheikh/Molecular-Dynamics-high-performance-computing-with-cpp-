cmake_minimum_required(VERSION 3.11)

project(molecular_dynamics)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)

set(MYPROJECT_HDRS
        # Add header files here
        ${CMAKE_CURRENT_SOURCE_DIR}/verlet.h
        types.h
        atoms.h lj_direct_summation.h)

set(MYPROJECT_SRCS
        # Add source files here
        ${CMAKE_CURRENT_SOURCE_DIR}/verlet.cpp
        tests/test_lj_direct_summation.cpp lj_direct_summation.cpp)

set(MYPROJECT_executable_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/verlet.cpp
        )


# Simulation code
add_executable(molecular_dynamics ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${MYPROJECT_executable_SRCS})
target_link_libraries(molecular_dynamics PUBLIC Eigen3::Eigen)
