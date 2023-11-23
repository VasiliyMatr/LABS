# Reference non-parallel program
add_executable(ref0 src/reference.cpp)
target_include_directories(ref0 PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(ref0 PRIVATE "${MPI_CXX_LIBRARIES}")

# Parallel program
add_executable(par0 src/parallel.cpp)
target_include_directories(par0 PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(par0 PRIVATE "${MPI_CXX_LIBRARIES}")
