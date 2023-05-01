#include <array>
#include <fstream>
#include <iostream>
#include <numbers>

#include <mpi.h>

#include <diff_solve.hpp>

#ifndef MPI_TRANS_EQ_HPP_INCL
#define MPI_TRANS_EQ_HPP_INCL

// Contains functionality for transfer equation MPI compute
namespace MPI_trans_eq {

using diff_compute::Float;

static constexpr MPI_Datatype get_float_datatype() noexcept {
  static_assert(std::is_same<Float, float>::value ||
                std::is_same<Float, double>::value);

  if (std::is_same<Float, float>::value) {
    return MPI_FLOAT;
  }

  return MPI_DOUBLE;
}

class Comm final {
  const int m_rank;
  const bool m_have_prev;
  const bool m_have_next;

public:
  Comm(int commsize, int rank)
      : m_rank(rank), m_have_prev(rank != 0),
        m_have_next((rank + 1) != commsize) {}

  Float get_prev_last() const noexcept {
    Float out = std::nan("");

    if (m_have_prev) {
      MPI_Recv(&out, 1, MPI_DOUBLE, m_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }

    return out;
  }
  Float get_next_first() const noexcept {
    Float out = std::nan("");

    if (m_have_next) {
      MPI_Recv(&out, 1, MPI_DOUBLE, m_rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }

    return out;
  }

  void set_curr_first(Float val) const noexcept {
    if (m_have_prev) {
      MPI_Send((void *)&val, 1, MPI_DOUBLE, m_rank - 1, 0, MPI_COMM_WORLD);
    }
  }
  void set_curr_last(Float val) const noexcept {
    if (m_have_next) {
      MPI_Send((void *)&val, 1, MPI_DOUBLE, m_rank + 1, 0, MPI_COMM_WORLD);
    }
  }
};

// Dump solutions grids for sequential x axis values ranges
inline void dump_grids(std::ofstream &out,
                       const std::vector<std::vector<Float>> &grids,
                       size_t t_size) {
  for (size_t t_idx = 0; t_idx != t_size; ++t_idx) {
    for (auto &&grid : grids) {
      size_t x_size = grid.size() / t_size;

      for (size_t x_idx = 0; x_idx != x_size; ++x_idx) {
        out << grid[t_idx * x_size + x_idx] << " ";
      }
    }
    out << std::endl;
  }
}

}; // namespace MPI_trans_eq

#endif // #ifndef MPI_TRANS_EQ_HPP_INCL
