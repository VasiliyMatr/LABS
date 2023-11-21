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

  init_a();

  // Workload
  // D = (-4, 2) => d = ('>', '<')
  auto start = MPI_Wtime();
  for (size_t i = 4; i != ISIZE; ++i) {
    for (size_t j = 0; j != JSIZE - 2; ++j) {
      a[i][j] = f(a[i - 4][j + 2]);
    }
  }
  auto end = MPI_Wtime();

  std::cout << "Workload time: " << end - start << std::endl;

  dump_a("ref.txt");

  MPI_Finalize();
  return 0;
}
