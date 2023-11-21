#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include <mpi.h>

#include "common.hpp"

auto get_job_range(int commsize, int rank) {
  bool is_last_job = rank == commsize - 1;

  size_t i_step = ISIZE / commsize;
  size_t i_begin = rank * i_step;
  size_t i_end = is_last_job ? ISIZE : i_begin + i_step;

  return std::make_pair(i_begin, i_end);
}

void calculate_a(int commsize, int rank) {
  auto [i_begin, i_end] = get_job_range(commsize, rank);

  for (size_t i = i_begin; i != i_end; ++i) {
    for (size_t j = 0; j != JSIZE; ++j) {
      a[i][j] = a_f(a[i][j]);
    }
  }

  if (rank != 0) {
    MPI_Send(a[i_begin], JSIZE * (i_end - i_begin), MPI_DOUBLE, 0, 0,
             MPI_COMM_WORLD);
  }
}

int main(int argc, char **argv) {
  int commsize = 0;
  int rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  init_ab();

  if (rank != 0) {
    calculate_a(commsize, rank);
    MPI_Finalize();
    return 0;
  }

  auto start = MPI_Wtime();

  std::vector<MPI_Request> requests{};
  for (int rank = 1; rank != commsize; ++rank) {
    auto [i_begin, i_end] = get_job_range(commsize, rank);

    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Irecv(a[i_begin], JSIZE * (i_end - i_begin), MPI_DOUBLE, rank, 0,
              MPI_COMM_WORLD, &request);

    requests.push_back(request);
  }

  // Workload a
  // D = (0, 0) => d = ('=', '=')
  calculate_a(commsize, rank);
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

  // Workload b
  // Not expensive => will calculate in one process
  // D = (1, -2) => d = ('>', '<')
  for (size_t i = 0; i != ISIZE - 1; ++i) {
    for (size_t j = 2; j != JSIZE; ++j) {
      b[i][j] = b_f(a[i + 1][j - 2]);
    }
  }
  auto end = MPI_Wtime();

  std::cout << end - start << std::endl;

  dump_b("par.txt");

  MPI_Finalize();
  return 0;
}
