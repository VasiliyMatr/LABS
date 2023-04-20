
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <numbers>

#include <diff_solve.hpp>

using diff_compute::Float;

class FuncT {
  const int m_rank;

  int m_row_idx = 1;

  Float func(Float t) const noexcept { return std::cos(t / 3); }

public:
  FuncT(int rank) : m_rank(rank) {}

  Float operator()(Float t) noexcept {
    if (m_rank == 0) {
      return func(t);
    }

    Float value = std::nan("");

    MPI_Recv(&value, 1, MPI_DOUBLE, m_rank - 1, m_row_idx++, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    return value;
  }
};

struct FuncX {
  Float operator()(Float x) const noexcept { return std::sin(x * x) + 1; }
};

struct FuncF {
  Float operator()(Float t, Float x) const noexcept {
    return t + std::exp(-x * x);
  }
};

class RowHandler {
  const int m_commsize;
  const int m_dest_rank;

  int row_idx = 1;

public:
  RowHandler(int commsize, int rank)
      : m_commsize(commsize), m_dest_rank(rank + 1) {}

  void operator()(Float value) noexcept {
    if (m_commsize == 1 || m_dest_rank == m_commsize) {
      return;
    }

    MPI_Send((void *)&value, 1, MPI_DOUBLE, m_dest_rank, row_idx++,
             MPI_COMM_WORLD);
  }
};

diff_compute::FloatRange
compute_rank_x_range(diff_compute::FloatRange total_range, int commsize,
                     int rank) {
  Float step = total_range.len() / commsize;

  return {total_range.min() + step * rank,
          total_range.min() + step * (rank + 1)};
}

using Solver = diff_compute::TransEqSolver<FuncT, FuncX, FuncF, RowHandler>;

int main(int argc, char **argv) {
  int commsize = 0;
  int my_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  diff_compute::FloatRange t_range(0., 5);
  diff_compute::FloatRange total_x_range(0., 5);
  Float dt = 0.005;
  Float dx = 0.005;

  diff_compute::FloatRange x_range =
      compute_rank_x_range(total_x_range, commsize, my_rank);

  // std::cout << "Rank " << my_rank << " computes for x = [" << x_range.min()
  //           << ", " << x_range.max() << "]" << std::endl;

  Solver::SolGridShape sol_grid_shape(t_range, x_range, dt, dx);

  Solver::TransEqFuncs funcs{FuncT{my_rank}, FuncX{}, FuncF{}};
  RowHandler row_handler(commsize, my_rank);

  auto sol = Solver(sol_grid_shape, funcs, row_handler).solve();

  if (my_rank != 0) {
    int sol_size = sol.size();

    // Send size and data
    MPI_Send((void *)&sol_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    MPI_Send((void *)sol.data(), sol.size(), MPI_DOUBLE, 0, 2,
             MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
  }

  std::vector<std::vector<Float>> sols = {std::move(sol)};

  for (int i = 1; i != commsize; ++i) {
    int sol_i_size = 0;

    MPI_Recv(&sol_i_size, 1, MPI_INT, i, 1, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    std::vector<Float> sol_i(sol_i_size, Float{0});

    MPI_Recv(sol_i.data(), sol_i_size, MPI_DOUBLE, i, 2,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    sols.emplace_back(std::move(sol_i));
  }

  size_t t_size = sol_grid_shape.t_size();

  for (size_t t_idx = 0; t_idx != t_size; ++t_idx) {
    for (auto &&sol : sols) {
      size_t x_size = sol.size() / t_size;
      for (size_t x_idx = 0; x_idx != x_size; ++x_idx) {
        std::cout << sol[t_idx * x_size + x_idx] << " ";
      }
      std::cout << std::endl;
    }
  }

  MPI_Finalize();
  return 0;
}
