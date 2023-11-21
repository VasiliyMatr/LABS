cmake_minimum_required(VERSION 3.11)

project(ParallelCycle VERSION 0.1 LANGUAGES CXX)

find_package(MPI REQUIRED)

# Compile commands export for language server
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -O2")

add_subdirectory(01)
add_subdirectory(03)
