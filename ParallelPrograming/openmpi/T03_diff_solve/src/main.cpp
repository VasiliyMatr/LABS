#include <cstdlib>
#include <functional>

#include <MPI_trans_eq.hpp>

using diff_compute::Float;

struct FuncT {
  Float operator()(Float t) const noexcept { return std::cos(t / 3); }
};

struct FuncX {
  Float operator()(Float x) const noexcept { return std::sin(x * x) + 1; }
};

Float func_f(Float t, Float x) noexcept {
  Float pow = std::pow(std::abs(x), std::cos(x) / 2);
  Float phase = std::pow(t * t + x * x, std::sin(x) / 2);

  return std::sin(phase) / std::exp(-pow);
}

struct FuncF {
  Float operator()(Float t, Float x) const noexcept {
    Float pow = std::pow(std::abs(x), std::cos(x) / 2);
    Float phase = std::pow(t * t + x * x, std::sin(x) / 2);

    return std::sin(phase) * std::exp(-pow);
  }
};

using Solver =
    diff_compute::TransEqSolver<FuncT, FuncX, FuncF, MPI_trans_eq::Comm>;

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
  diff_compute::FloatRange x_range(0., 5);

  bool need_dump = argc == 4;

  assert_with_msg(argc == 3 || need_dump, "Invalid arguments");

  Float dt = std::atof(argv[1]);
  Float dx = std::atof(argv[2]);

  assert_with_msg(dt != 0 && dx != 0, "Invalid arguments");

  Solver::GridShape grid_shape(t_range, x_range, dt, dx);
  Solver::GridShape sub_grid_shape =
      grid_shape.get_x_subgrid(commsize, my_rank);

  std::cout << "Rank " << my_rank << " computes for x = ["
            << sub_grid_shape.x_range().min() << ", "
            << sub_grid_shape.x_range().max() << "]" << std::endl;

  double t_begin = MPI_Wtime();

  FuncT t{};
  FuncX x{};
  FuncF f{};

  MPI_trans_eq::Comm comm{commsize, my_rank};

  auto grid = Solver(t, x, f, sub_grid_shape, comm).solve();

  if (my_rank != 0) {
    size_t grid_size = grid.size();

    // Send grid size and data
    MPI_Send((void *)&grid_size, 1, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
    MPI_Send((void *)grid.data(), grid_size, MPI_trans_eq::get_float_datatype(),
             0, 2, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
  }

  double t_no_send = MPI_Wtime();
  std::cout << "No send compute time: " << t_no_send - t_begin << "sec"
            << std::endl;

  std::vector<std::vector<Float>> grids{};
  grids.emplace_back(std::move(grid));

  for (int i = 1; i != commsize; ++i) {
    size_t grid_i_size = 0;

    MPI_Recv(&grid_i_size, 1, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    std::vector<Float> grid_i(grid_i_size, Float{0});

    MPI_Recv(grid_i.data(), grid_i_size, MPI_trans_eq::get_float_datatype(), i,
             2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    grids.emplace_back(std::move(grid_i));
  }

  double t_end = MPI_Wtime();

  std::cout << "Total compute time: " << t_end - t_begin << "sec" << std::endl;

  if (need_dump) {
    std::ofstream out;
    out.open(argv[3]);
    assert_with_msg(out.is_open(), "Can not open output file");
    MPI_trans_eq::dump_grids(out, grids, grid_shape.t_size());
  }

  MPI_Finalize();
  return 0;
}
