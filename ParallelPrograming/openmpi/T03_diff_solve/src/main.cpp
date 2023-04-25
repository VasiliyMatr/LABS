#include <cstdlib>

#include <MPI_trans_eq.hpp>

using diff_compute::Float;

struct FuncT {
  Float operator()(Float t) const noexcept { return std::cos(t / 3); }
};

struct FuncX {
  Float operator()(Float x) const noexcept { return std::sin(x * x) + 1; }
};

struct FuncF {
  Float operator()(Float t, Float x) const noexcept {
    return t + std::exp(-x * x);
  }
};

static constexpr size_t BUFF_SIZE = 32;

using Solver = MPI_trans_eq::TransEqSolver<FuncT, FuncX, FuncF, BUFF_SIZE>;
using MPI_FuncT = MPI_trans_eq::FuncT<FuncT, BUFF_SIZE>;
using MPI_RowHandler = MPI_trans_eq::RowHandler<BUFF_SIZE>;

void assert_with_msg(bool cond, const char *msg, int exit_code = -1) {
  if (!cond) {
    std::cout << msg << std::endl;
    exit(exit_code);
  }
}

int main(int argc, char **argv) {
  int commsize = 0;
  int my_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  diff_compute::FloatRange t_range(0., 5);
  diff_compute::FloatRange total_x_range(0., 5);

  bool need_dump = argc == 4;

  assert_with_msg(argc == 3 || need_dump, "Invalid arguments");

  Float dt = std::atof(argv[1]);
  Float dx = std::atof(argv[2]);

  assert_with_msg(dt != 0 && dx != 0, "Invalid arguments");

  auto x_range =
      MPI_trans_eq::compute_rank_x_range(total_x_range, commsize, my_rank);

  std::cout << "Rank " << my_rank << " computes for x = [" << x_range.min()
            << ", " << x_range.max() << "]" << std::endl;

  Solver::SolGridShape sol_grid_shape(t_range, x_range, dt, dx);

  Solver::TransEqFuncs funcs{MPI_FuncT{my_rank}, FuncX{}, FuncF{}};
  MPI_RowHandler row_handler{commsize, my_rank, sol_grid_shape.t_size()};

  auto grid = Solver(sol_grid_shape, funcs, row_handler).solve();

  if (my_rank != 0) {
    int grid_size = grid.size();

    // Send grid size and data
    MPI_Send((void *)&grid_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    MPI_Send((void *)grid.data(), grid_size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
  }

  std::vector<std::vector<Float>> grids = {std::move(grid)};

  for (int i = 1; i != commsize; ++i) {
    int grid_i_size = 0;

    MPI_Recv(&grid_i_size, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::vector<Float> grid_i(grid_i_size, Float{0});

    MPI_Recv(grid_i.data(), grid_i_size, MPI_DOUBLE, i, 2, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    grids.emplace_back(std::move(grid_i));
  }

  if (need_dump) {
    std::ofstream out;
    out.open(argv[3]);
    assert_with_msg(out.is_open(), "Can not open output file");
    MPI_trans_eq::dump_grids(out, grids, sol_grid_shape.t_size());
  }

  MPI_Finalize();
  return 0;
}
