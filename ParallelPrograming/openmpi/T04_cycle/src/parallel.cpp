#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <utility>
#include <vector>

#include <mpi.h>

#include "common.hpp"

static constexpr size_t IDIST = 4;
static constexpr size_t JDIST = 2;

// Calculate initial job bounds
auto calc_job_bounds(int commsize, int rank) {
  bool is_last_job = rank == commsize - 1;

  size_t job_step = JSIZE / commsize;
  size_t job_begin = rank * job_step;
  size_t job_end = is_last_job ? JSIZE : job_begin + job_step;

  return std::make_pair(job_begin, job_end);
}

// Collect i-th result of rank-th job
void collect(int commsize, int rank, size_t i,
             const std::vector<double> &data) {
  auto [job_begin, job_end] = calc_job_bounds(commsize, rank);
  size_t job_size = job_end - job_begin;

  size_t offset = (i / IDIST) * JDIST;

  for (size_t data_i = 0; data_i != IDIST; ++data_i) {
    for (size_t j = 0; j != job_size; ++j) {
      size_t shifted_j = (job_begin + j - offset + JSIZE) % JSIZE;

      a[i + data_i][shifted_j] = data[data_i * job_size + j];
    }
  }
}

// Compute my job
void compute(int commsize, int my_rank) {
  auto [job_begin, job_end] = calc_job_bounds(commsize, my_rank);
  size_t job_size = job_end - job_begin;

  std::vector<double> data{};
  data.reserve(job_size * IDIST);

  for (size_t i = 0; i != IDIST; ++i) {
    for (size_t j = 0; j != job_size; ++j) {
      data.push_back(a[i][j + job_begin]);
    }
  }

  for (size_t i = IDIST, offset = JDIST; i != ISIZE; ++i) {
    size_t data_i = i % IDIST;

    for (size_t j = 0; j != job_size; ++j) {
      size_t shifted_j = (job_begin + j - offset + JSIZE) % JSIZE;

      if (shifted_j >= JSIZE - JDIST) {
        data[data_i * job_size + j] = a[i][shifted_j];
      } else {
        data[data_i * job_size + j] = f(data[data_i * job_size + j]);
      }
    }

    if (data_i == IDIST - 1) {
      offset += JDIST;
      size_t tag = i - data_i;

      if (my_rank == 0) {
        collect(commsize, 0, tag, data);
      } else {
        MPI_Send(data.data(), data.size(), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
      }
    }
  }
}

// Collect all results
void collect(int commsize) {
  auto [last_job_begin, last_job_end] = calc_job_bounds(commsize, commsize - 1);
  size_t max_job_size = last_job_end - last_job_begin;
  std::vector<double> data(max_job_size * IDIST, std::nan(""));

  for (size_t i = IDIST; i < ISIZE; i += 4) {
    for (int rank = 1; rank != commsize; ++rank) {
      auto [job_begin, job_end] = calc_job_bounds(commsize, rank);
      size_t job_size = job_end - job_begin;

      MPI_Recv(data.data(), job_size * IDIST, MPI_DOUBLE, rank, i,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      collect(commsize, rank, i, data);
    }
  }
}

int main(int argc, char **argv) {
  int commsize = 0;
  int my_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  init_a();

  // Workload
  // D = (-4, 2) => d = ('>', '<')
  if (my_rank != 0) {
    compute(commsize, my_rank);

    MPI_Finalize();
    return 0;
  }

  auto start = MPI_Wtime();
  compute(commsize, my_rank);
  collect(commsize);
  auto end = MPI_Wtime();

  std::cout << "Workload time: " << end - start << std::endl;

  dump_a("par.txt");

  MPI_Finalize();
  return 0;
}
