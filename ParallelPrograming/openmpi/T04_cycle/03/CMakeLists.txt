# Reference non-parallel program
add_executable(ref3 src/reference.cpp)
target_include_directories(ref3 PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(ref3 PRIVATE "${MPI_CXX_LIBRARIES}")

# Parallel program
add_executable(par3 src/parallel.cpp)
target_include_directories(par3 PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(par3 PRIVATE "${MPI_CXX_LIBRARIES}")
