#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "common.hpp"

int main(int argc, char **argv) {
  int commsize = 0;
  int my_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  assert(commsize == 1);

  init_ab();

  auto start = MPI_Wtime();
  // Workload a
  // D = (0, 0) => d = ('=', '=')
  for (size_t i = 0; i != ISIZE; ++i) {
    for (size_t j = 0; j != JSIZE; ++j) {
      a[i][j] = a_f(a[i][j]);
    }
  }
  auto a_end = MPI_Wtime();

  // Workload b
  // D = (1, -2) => d = ('>', '<')
  for (size_t i = 0; i != ISIZE - 1; ++i) {
    for (size_t j = 2; j != JSIZE; ++j) {
      b[i][j] = b_f(a[i + 1][j - 2]);
    }
  }
  auto end = MPI_Wtime();

  std::cout << "A time: " << a_end - start << std::endl;
  std::cout << "Total time: " << end - start << std::endl;

  dump_b("ref.txt");

  MPI_Finalize();
  return 0;
}
