cmake_minimum_required(VERSION 3.11)

project(ParallelCycle VERSION 0.1 LANGUAGES CXX)

find_package(MPI REQUIRED)

# Compile commands export for language server
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -O2")

# Reference non-parallel program
add_executable(reference src/reference.cpp)
target_include_directories(reference PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(reference PRIVATE "${MPI_CXX_LIBRARIES}")

# Parallel program
add_executable(parallel src/parallel.cpp)
target_include_directories(parallel PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(parallel PRIVATE "${MPI_CXX_LIBRARIES}")
