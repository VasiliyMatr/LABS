
#include <cassert>
#include <cmath>
#include <iterator>
#include <vector>

#include <iostream>

namespace diff_compute {

using Float = double;

class FloatRange final {
  Float m_min = std::nan("");
  Float m_max = std::nan("");
  Float m_len = std::nan("");

public:
  FloatRange() {}

  FloatRange(Float min, Float max) : m_min(min), m_max(max), m_len(max - min) {
    assert(min < max && "Invalid range bounds");
  }

  Float min() const noexcept { return m_min; }
  Float max() const noexcept { return m_max; }
  Float len() const noexcept { return m_len; }
};

template <class FuncT, class FuncX, class FuncF, class RowHandler>
struct TransEqSolver final {
  class SolGridShape final {
    FloatRange m_t_range = FloatRange();
    FloatRange m_x_range = FloatRange();

    Float m_dt = std::nan("");
    Float m_dx = std::nan("");

    size_t m_t_size = 0;
    size_t m_x_size = 0;

  public:
    SolGridShape(FloatRange t_range, FloatRange x_range, Float dt, Float dx)
        : m_t_range(t_range), m_x_range(x_range), m_dt(dt), m_dx(dx),
          m_t_size(t_range.len() / dt), m_x_size(x_range.len() / dx) {
      assert(dx >= dt && "The solution may diverge for dx >= dt");
      assert(m_t_size != 0 && m_x_size != 0 && "Invalid SolGrid shape");
    }

    // Get t axis FloatRange
    FloatRange t_range() const noexcept { return m_t_range; }
    // Get x axis FloatRange
    FloatRange x_range() const noexcept { return m_x_range; }

    // Get t axis split step
    Float dt() const noexcept { return m_dt; }
    // Get x axis split step
    Float dx() const noexcept { return m_dx; }

    // Get t axis elements number
    size_t t_size() const noexcept { return m_t_size; }
    // Get x axis elements number
    size_t x_size() const noexcept { return m_x_size; }

    // Calculate t coordinate for given t index
    Float t_coord(size_t t_idx) const noexcept {
      Float t_delta = m_dt * t_idx;
      return m_t_range.min() + t_delta;
    }
    // Calculate x coordinate for given x index
    Float x_coord(size_t x_idx) const noexcept {
      Float x_delta = m_dx * x_idx;
      return m_x_range.min() + x_delta;
    }
  };

  struct SolGridPoint final {
    Float t = std::nan("");
    Float x = std::nan("");
    Float U = std::nan("");
  };

  class TransEqFuncs final {
    FuncT m_t;
    FuncX m_x;
    FuncF m_f;

  public:
    TransEqFuncs(FuncT t, FuncX x, FuncF f) : m_t(t), m_x(x), m_f(f) {}

    FuncT &t() noexcept { return m_t; }
    FuncX &x() noexcept { return m_x; }
    FuncF &f() noexcept { return m_f; }
  };

private:
  template <class UIt> class RowIt final {
    UIt m_U_it{};
    const SolGridShape &m_sol_grid_shape;
    const size_t m_t_idx;
    size_t m_x_idx = 0;

  public:
    RowIt(UIt U_begin, const SolGridShape &sol_grid_shape, size_t t_idx)
        : m_U_it(U_begin + t_idx * sol_grid_shape.x_size()),
          m_sol_grid_shape(sol_grid_shape), m_t_idx(t_idx) {}

    SolGridPoint operator*() const noexcept {
      Float t = m_sol_grid_shape.t_coord(m_t_idx);
      Float x = m_sol_grid_shape.x_coord(m_x_idx);
      return {t, x, *m_U_it};
    }
    RowIt &operator++() noexcept {
      ++m_U_it;
      ++m_x_idx;
      return *this;
    }
    bool operator==(const RowIt &other) const noexcept {
      return m_U_it == other.m_U_it;
    }
  };

  SolGridShape m_sol_grid_shape;
  TransEqFuncs m_trans_eq_funcs;
  RowHandler m_row_handler;

public:
  TransEqSolver(const SolGridShape &sol_grid_shape,
                const TransEqFuncs &trans_eq_funcs,
                const RowHandler &row_handler)
      : m_sol_grid_shape(sol_grid_shape), m_trans_eq_funcs(trans_eq_funcs),
        m_row_handler(row_handler) {}

private:
  Float angle_method_step(Float U1, Float U2, Float Fval) const noexcept {
    return U2 +
           m_sol_grid_shape.dt() * (Fval + (U1 - U2) / m_sol_grid_shape.dx());
  }

  template <class PointIt, class FloatIt>
  void compute_row_with_angle_method(PointIt prev_row_it, PointIt prev_row_end,
                                     FloatIt out_it) {
    SolGridPoint P1 = *prev_row_it;
    ++prev_row_it;
    SolGridPoint P2{};

    for (; prev_row_it != prev_row_end; ++prev_row_it, ++out_it) {
      P2 = *prev_row_it;
      *out_it = angle_method_step(P1.U, P2.U, m_trans_eq_funcs.f()(P2.t, P2.x));

      P1 = P2;
    }
  }

  Float cross_method_step(Float U1, Float U2, Float U3,
                          Float Fval) const noexcept {
    return U2 + m_sol_grid_shape.dt() *
                    (Fval * 2 + (U1 - U3) / m_sol_grid_shape.dx());
  }

  template <class PointIt, class FloatIt>
  void compute_row_with_cross_method(PointIt row1_it, PointIt row1_end,
                                     PointIt row2_it, FloatIt out_it) {
    SolGridPoint P1 = *row1_it;
    ++row1_it;
    SolGridPoint next_P1 = *row1_it;
    ++row1_it;
    SolGridPoint P2{};
    SolGridPoint P3{};

    for (; row1_it != row1_end; ++row1_it, ++row2_it) {
      P2 = *row2_it;
      P3 = *row1_it;

      *out_it =
          cross_method_step(P1.U, P2.U, P3.U, m_trans_eq_funcs.f()(P1.t, P2.x));

      P1 = next_P1;
      next_P1 = P3;
      ++out_it;
    }

    *out_it = angle_method_step(P1.U, next_P1.U,
                                m_trans_eq_funcs.f()(next_P1.t, next_P1.x));
  }

public:
  std::vector<Float> solve() {
    std::vector<Float> grid_data{};
    size_t total_grid_size =
        m_sol_grid_shape.t_size() * m_sol_grid_shape.x_size();
    grid_data.reserve(total_grid_size);

    size_t x_size = m_sol_grid_shape.x_size();
    size_t t_size = m_sol_grid_shape.t_size();

    // Compute first row - with given border function
    for (size_t x_idx = 0; x_idx != x_size; ++x_idx) {
      Float x_coord = m_sol_grid_shape.x_coord(x_idx);
      grid_data.push_back(m_trans_eq_funcs.x()(x_coord));
    }

    // Compute other rows - with angular method
    for (size_t t_idx = 1; t_idx != t_size; ++t_idx) {
      Float t_coord = m_sol_grid_shape.t_coord(t_idx);

      // Compute first value with given border function
      Float border_val = m_trans_eq_funcs.t()(t_coord);
      grid_data.push_back(border_val);

      // Prepare iterators
      RowIt prev_row1_it(grid_data.begin(), m_sol_grid_shape, t_idx - 1);
      RowIt prev_row1_end(grid_data.begin(), m_sol_grid_shape, t_idx);
      auto curr_row_it = std::inserter(grid_data, grid_data.end());

      if (t_idx == 1) {
        compute_row_with_angle_method(prev_row1_it, prev_row1_end, curr_row_it);
      } else {
        RowIt prev_row2_it(grid_data.begin(), m_sol_grid_shape, t_idx - 2);
        ++prev_row2_it;

        compute_row_with_cross_method(prev_row1_it, prev_row1_end, prev_row2_it,
                                      curr_row_it);
      }

      // Pass row last value to given handler
      m_row_handler(grid_data.back());
    }

    return grid_data;
  }
};

} // namespace diff_compute
