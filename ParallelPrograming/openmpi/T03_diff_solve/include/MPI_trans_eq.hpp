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

// Time boundary function for transfer equation that can
// calculate values with given FuncT or recieve those values
// from MPI functions
// - Buffering optimization is implemented
template <class Func, size_t buff_size> class FuncT {
  static_assert(buff_size != 0);

  const int m_rank;
  const Func m_func{};

  std::array<Float, buff_size> m_buff{};
  int m_request_idx = 0;

public:
  FuncT(int rank) : m_rank(rank) {}

  // Uses provided function for executor with zero rank
  // Otherwise recieves values form previous rank via MPI functions
  Float operator()(Float t) noexcept {
    if (m_rank == 0) {
      return m_func(t);
    }

    if (m_request_idx % buff_size == 0) {
      MPI_Recv(m_buff.data(), buff_size, MPI_DOUBLE, m_rank - 1, MPI_ANY_TAG,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return m_buff[m_request_idx++ % buff_size];
  }
};

// Transfer equation last value in rows handler
// Sends last value in rows to the next rank executor
// - Buffering optimization is implemented
template <size_t buff_size> class RowHandler {
  static_assert(buff_size != 0);

  const int m_commsize;
  const int m_dest_rank;
  const size_t m_t_size;

  std::array<Float, buff_size> m_buff{};
  size_t m_request_idx = 0;

  // Fill m_buff with NaN values
  void reset_buff() noexcept {
    for (size_t i = 0; i != buff_size; ++i) {
      m_buff[i] = std::nan("");
    }
  }

public:
  RowHandler(int commsize, int rank, size_t t_size)
      : m_commsize(commsize), m_dest_rank(rank + 1), m_t_size(t_size) {
    reset_buff();
  }

  void operator()(Float value) noexcept {
    if (m_dest_rank == m_commsize) {
      return;
    }

    m_buff[m_request_idx % buff_size] = value;

    bool buff_full = (m_request_idx + 1) % buff_size == 0;
    bool last_request = (m_request_idx + 1) == m_t_size - 1;

    if (buff_full || last_request) {
      MPI_Send((void *)m_buff.data(), buff_size, MPI_DOUBLE, m_dest_rank, 0,
               MPI_COMM_WORLD);
      reset_buff();
    }

    ++m_request_idx;
  }
};

template <class T, class X, class F, size_t buff_size>
using TransEqSolver = diff_compute::TransEqSolver<FuncT<T, buff_size>, X, F,
                                                  RowHandler<buff_size>>;

// Calculate x axis range to be processed by given rank
inline diff_compute::FloatRange
compute_rank_x_range(diff_compute::FloatRange total_range, int commsize,
                     int rank) {
  Float step = total_range.len() / commsize;

  return {total_range.min() + step * rank,
          total_range.min() + step * (rank + 1)};
}

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
      out << std::endl;
    }
  }
}

}; // namespace MPI_trans_eq

#endif // #ifndef MPI_TRANS_EQ_HPP_INCL
