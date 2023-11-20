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

  size_t job_step = JSIZE / (commsize - 1);
  size_t job_begin = (rank - 1) * job_step;
  size_t job_end = is_last_job ? JSIZE : job_begin + job_step;

  return std::make_pair(job_begin, job_end);
}

// Collect all results
void collect(int commsize) {
  std::vector<MPI_Request> requests(commsize - 1, MPI_REQUEST_NULL);

  int copy_end = -1;

  for (size_t i = IDIST, offset = JDIST; i < ISIZE; ++i) {
    for (int rank = 1; rank != commsize; ++rank) {
      auto [job_begin, job_end] = calc_job_bounds(commsize, rank);
      size_t job_size = job_end - job_begin;

      size_t shifted_begin = (job_begin - offset + JSIZE) % JSIZE;

      MPI_Irecv(a[i] + shifted_begin, job_size, MPI_DOUBLE, rank, i,
                MPI_COMM_WORLD, &requests[rank - 1]);

      if (shifted_begin + job_size > JSIZE) {
        copy_end = shifted_begin + job_size - JSIZE;
      }
    }

    for (auto &&request : requests) {
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    for (int j = 0; j < copy_end; ++j) {
      a[i][j] = a[i + 1][j];
    }

    copy_end = -1;
    if (i % IDIST == IDIST - 1) {
      offset += JDIST;
    }
  }
}

// Compute my job
void compute(int commsize, int my_rank) {
  auto [job_begin, job_end] = calc_job_bounds(commsize, my_rank);
  size_t job_size = job_end - job_begin;

  std::array<std::vector<double>, ISIZE> data{};
  std::array<MPI_Request, IDIST> requests{};

  // Init start data
  for (auto &&idata : data) {
    idata.reserve(job_size);
  }

  for (size_t i = 0; i != IDIST; ++i) {
    for (size_t j = job_begin; j != job_end; ++j) {
      data[i].push_back(a[i][j]);
    }

    requests[i] = MPI_REQUEST_NULL;
  }

  // Compute and send
  for (size_t i = IDIST, offset = JDIST; i != ISIZE; ++i) {
    size_t imod = i % IDIST;
    auto &&idata = data[imod];

    for (size_t j = 0; j != job_size; ++j) {
      size_t shifted_j = (job_begin + j + JSIZE - offset) % JSIZE;

      if (shifted_j >= JSIZE - JDIST) {
        idata[j] = a[i][shifted_j];
      } else {
        idata[j] = f(idata[j]);
      }
    }

    MPI_Isend(idata.data(), idata.size(), MPI_DOUBLE, 0, i, MPI_COMM_WORLD,
              &requests[imod]);

    if (imod == IDIST - 1) {
      // Increment offset
      offset += JDIST;

      // Wait for all requests
      for (auto &&request : requests) {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }
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
  collect(commsize);
  auto end = MPI_Wtime();

  std::cout << "Workload time: " << end - start << std::endl;

  dump_a("par.txt");

  MPI_Finalize();
  return 0;
}
