cmake_minimum_required(VERSION 3.30)
project(hw1_comp_geo_problem_1)

set(CMAKE_CXX_STANDARD 20)

# Find CGAL package
find_package(CGAL REQUIRED)

# Include CGAL headers and link against CGAL libraries
include(${CGAL_USE_FILE})

set(SOURCE_FILES main.cpp)

# Add executables
add_executable(hw1_comp_geo_problem_1 main.cpp
        convexImplementations.h
        polygonalTriangulationImplementations.h
        kdtree.h
        rangetree.h)

# Link against CGAL libraries
target_link_libraries(hw1_comp_geo_problem_1 ${CGAL_LIBS})
