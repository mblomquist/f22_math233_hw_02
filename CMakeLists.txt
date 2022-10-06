cmake_minimum_required(VERSION 3.22)
project(f22_math233_hw_02)

set(CMAKE_CXX_STANDARD 14)

add_executable(f22_math233_hw_02 src/main.cpp src/grid/Grid2d.cpp src/grid/Grid2d.h src/advection/SL_method.cpp src/advection/SL_method.h src/tools/math_tools.cpp src/tools/math_tools.h)
