#include <array>
#include <cassert>
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

class JobConfig final {
  size_t m_commsize = 0;
  size_t m_rank = 0;

  bool m_is_last_job = m_rank == m_commsize - 1;

  size_t m_j_step = (JSIZE - JDIST) / m_commsize;
  size_t m_j_begin = m_rank * m_j_step;
  size_t m_j_end = m_is_last_job ? JSIZE - JDIST : m_j_begin + m_j_step;

public:
  JobConfig(int commsize, int rank) : m_commsize(commsize), m_rank(rank) {
    assert(commsize > 0);
    assert(rank >= 0);
    assert(commsize > rank);
    assert(m_j_step >= JDIST * 2);
  }

  auto commsize() const { return m_commsize; }
  auto rank() const { return m_rank; }
  auto is_last_job() const { return m_is_last_job; }
  auto j_begin() const { return m_j_begin; }
  auto j_end() const { return m_j_end; }
  auto j_size() const { return m_j_end - m_j_begin; }
};

int dep_tag(size_t i) { return i; }

int collect_tag(size_t i) { return i + ISIZE; }

void calculate(JobConfig &cfg) {
  std::vector<MPI_Request> collect_requests{};
  if (cfg.rank() != 0) {
    collect_requests.reserve(ISIZE - IDIST);
  }

  size_t job_deps_num = ISIZE - IDIST * 2;

  std::vector<MPI_Request> prev_deps_requests{};
  if (cfg.rank() != 0) {
    prev_deps_requests.reserve(job_deps_num);
  }

  std::vector<MPI_Request> job_deps_requests{};
  if (!cfg.is_last_job()) {
    job_deps_requests.reserve(job_deps_num);

    for (size_t i = IDIST * 2; i != ISIZE; ++i) {
      MPI_Request r = MPI_REQUEST_NULL;
      MPI_Irecv(a[i - IDIST] + cfg.j_end(), JDIST, MPI_DOUBLE, cfg.rank() + 1,
                dep_tag(i), MPI_COMM_WORLD, &r);
      job_deps_requests.push_back(r);
    }
  }

  for (size_t i = IDIST; i != ISIZE; ++i) {
    auto j = cfg.j_begin();

    for (size_t end = j + JDIST; j != end; ++j) {
      a[i][j] = f(a[i - IDIST][j + JDIST]);
    }

    if (cfg.rank() != 0 && i < ISIZE - IDIST) {
      MPI_Request r = MPI_REQUEST_NULL;
      MPI_Isend(a[i] + cfg.j_begin(), JDIST, MPI_DOUBLE, cfg.rank() - 1,
                dep_tag(i + IDIST), MPI_COMM_WORLD, &r);

      prev_deps_requests.push_back(r);
    }

    // if (!cfg.is_last_job() && i >= IDIST * 2) {
    //   MPI_Wait(&job_deps_requests[i - IDIST * 2], MPI_STATUS_IGNORE);
    // }

    for (size_t end = cfg.j_end(); j != end; ++j) {
      a[i][j] = f(a[i - IDIST][j + JDIST]);
    }

    if (cfg.rank() != 0) {
      MPI_Request r = MPI_REQUEST_NULL;
      MPI_Isend(a[i] + cfg.j_begin(), cfg.j_size(), MPI_DOUBLE, 0,
                collect_tag(i), MPI_COMM_WORLD, &r);
      collect_requests.push_back(r);
    }

    if (!cfg.is_last_job() && i % IDIST == IDIST - 1 && i != ISIZE - 1) {
      MPI_Waitall(IDIST, &job_deps_requests[i + 1 - IDIST * 2],
                  MPI_STATUS_IGNORE);
    }
  }

  MPI_Waitall(collect_requests.size(), collect_requests.data(),
              MPI_STATUS_IGNORE);
  MPI_Waitall(prev_deps_requests.size(), prev_deps_requests.data(),
              MPI_STATUS_IGNORE);
}

std::vector<MPI_Request> collect_jobs(JobConfig &cfg) {
  auto commsize = cfg.commsize();

  assert(cfg.rank() == 0);

  std::vector<MPI_Request> collect_requests{};
  collect_requests.reserve((ISIZE - IDIST) * (commsize - 1));

  for (size_t rank = 1; rank != commsize; ++rank) {
    JobConfig rank_cfg(commsize, rank);

    for (size_t i = IDIST; i != ISIZE; ++i) {
      MPI_Request r = MPI_REQUEST_NULL;
      MPI_Irecv(a[i] + rank_cfg.j_begin(), rank_cfg.j_size(), MPI_DOUBLE, rank,
                collect_tag(i), MPI_COMM_WORLD, &r);
      collect_requests.push_back(r);
    }
  }

  return collect_requests;
}

int main(int argc, char **argv) {
  int commsize = 0;
  int my_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  init_a();

  JobConfig cfg(commsize, my_rank);

  // Workload
  // D = (-4, 2) => d = ('>', '<')
  if (my_rank != 0) {
    calculate(cfg);

    MPI_Finalize();
    return 0;
  }

  auto start = MPI_Wtime();
  auto collect_requests = collect_jobs(cfg);
  calculate(cfg);

  MPI_Waitall(collect_requests.size(), collect_requests.data(),
              MPI_STATUS_IGNORE);
  auto end = MPI_Wtime();

  std::cout << "Workload time: " << end - start << std::endl;

  dump_a("par.txt");

  MPI_Finalize();
  return 0;
}
