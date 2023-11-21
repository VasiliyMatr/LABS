# Reference non-parallel program
add_executable(ref1 src/reference.cpp)
target_include_directories(ref1 PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(ref1 PRIVATE "${MPI_CXX_LIBRARIES}")

# Parallel program
add_executable(par1 src/parallel.cpp)
target_include_directories(par1 PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(par1 PRIVATE "${MPI_CXX_LIBRARIES}")
