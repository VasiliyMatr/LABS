
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <vector>

#include <iostream>

#ifndef DIFF_SOLVE_HPP_INCL
#define DIFF_SOLVE_HPP_INCL

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

// FuncF:
// - double FuncF::operator()(double x, double t) const
// FuncX:
// - double FuncX::operator()(double x) const
// FuncT:
// - double FuncT::operator()(double t) const
// Comm:
// - double Comm::get_prev_last()
// - double Comm::get_next_first()
// - double Comm::set_curr_first()
// - double Comm::set_curr_last()
template <class FuncT, class FuncX, class FuncF, class Comm>
struct TransEqSolver final {
  class GridShape final {
    FloatRange m_t_range = FloatRange();
    FloatRange m_x_range = FloatRange();

    Float m_dt = std::nan("");
    Float m_dx = std::nan("");

    size_t m_t_size = 0;
    size_t m_x_size = 0;

  public:
    GridShape(FloatRange t_range, FloatRange x_range, Float dt, Float dx)
        : m_t_range(t_range), m_x_range(x_range), m_dt(dt), m_dx(dx),
          m_t_size(1 + t_range.len() / dt), m_x_size(1 + x_range.len() / dx) {
      assert(dx >= dt && "The solution may diverge for dx < dt");
      assert(m_t_size != 0 && m_x_size != 0 && "Invalid GridShape");
    }

    GridShape get_x_subgrid(size_t subgrids_num,
                            size_t subgrid_idx) const noexcept {
      size_t x_size_step = m_x_size / subgrids_num;
      size_t subgrid_x_size = x_size_step;
      if (subgrids_num == subgrid_idx + 1) {
        subgrid_x_size += m_x_size % subgrids_num;
      }

      size_t subgrid_x_min_idx = x_size_step * subgrid_idx;
      size_t subgrid_x_max_idx = subgrid_x_min_idx + subgrid_x_size - 1;

      return GridShape(m_t_range,
                       {subgrid_x_min_idx * m_dx, subgrid_x_max_idx * m_dx},
                       m_dt, m_dx);
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
    Float t_coord(long long t_idx) const noexcept {
      Float t_delta = m_dt * t_idx;
      return m_t_range.min() + t_delta;
    }
    // Calculate x coordinate for given x index
    Float x_coord(long long x_idx) const noexcept {
      Float x_delta = m_dx * x_idx;
      return m_x_range.min() + x_delta;
    }
  };

private:
  // Represents 3D point: (t, x, U(x, t))
  struct Point final {
    Float t = std::nan("");
    Float x = std::nan("");
    Float u = std::nan("");

    bool is_valid() const noexcept {
      return !std::isnan(t) && !std::isnan(x) && !std::isnan(u);
    }
  };

  // Angle method input points
  using AnglePoints = std::array<Point, 2>;
  // Cross method input points
  using CrossPoints = std::array<Point, 3>;

  // @brief
  // Iterates trough grid row points
  template <std::input_iterator ValueIt> struct RowIt final {
    using difference_type = std::ptrdiff_t;
    using value_type = Point;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::input_iterator_tag;

    ValueIt m_value_it{};

    const GridShape m_grid_shape;
    const size_t m_t_idx;
    size_t m_x_idx = 0;

  public:
    RowIt(const ValueIt &grid_begin, const GridShape &grid_shape, size_t t_idx)
        : m_value_it(grid_begin + t_idx * grid_shape.x_size()),
          m_grid_shape(grid_shape), m_t_idx(t_idx) {}

    const GridShape &grid_shape() const noexcept { return m_grid_shape; }
    Float t_coord() const noexcept { return m_grid_shape.t_coord(m_t_idx); }

    value_type operator*() const noexcept {
      Float t = m_grid_shape.t_coord(m_t_idx);
      Float x = m_grid_shape.x_coord(m_x_idx);
      return {t, x, *m_value_it};
    }
    RowIt &operator++() noexcept {
      ++m_value_it;
      ++m_x_idx;
      return *this;
    }
    bool operator==(const RowIt &other) const noexcept {
      return m_value_it == other.m_value_it;
    }
  };

  template <std::input_iterator ValueIt> struct AngleIt final {
    using difference_type = std::ptrdiff_t;
    using value_type = AnglePoints;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::input_iterator_tag;

  private:
    RowIt<ValueIt> m_prev_row_it;
    Point m_prev_point;

  public:
    AngleIt(const RowIt<ValueIt> &prev_row_it) : m_prev_row_it(prev_row_it) {}

    value_type operator*() noexcept {
      if (!m_prev_point.is_valid()) {
        m_prev_point = *m_prev_row_it;
        ++m_prev_row_it;
      }

      return {m_prev_point, *m_prev_row_it};
    }
    AngleIt &operator++() noexcept {
      m_prev_point = *m_prev_row_it;
      ++m_prev_row_it;
      return *this;
    }
    bool operator==(const AngleIt &other) const noexcept {
      return m_prev_row_it == other.m_prev_row_it;
    }
  };

  template <std::input_iterator ValueIt> struct CrossIt final {
    using difference_type = std::ptrdiff_t;
    using value_type = CrossPoints;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::input_iterator_tag;

  private:
    RowIt<ValueIt> m_prev_row1_it;
    RowIt<ValueIt> m_prev_row2_it;

    Point m_p1{};
    Point m_next_p1{};

  public:
    CrossIt(const RowIt<ValueIt> &prev_row1_it,
            const RowIt<ValueIt> &prev_row2_it)
        : m_prev_row1_it(prev_row1_it), m_prev_row2_it(prev_row2_it) {
      ++m_prev_row2_it;
    }

    value_type operator*() noexcept {
      if (!m_p1.is_valid()) {
        m_p1 = *m_prev_row1_it;
        ++m_prev_row1_it;
        m_next_p1 = *m_prev_row1_it;
        ++m_prev_row1_it;
      }

      return {m_p1, *m_prev_row2_it, *m_prev_row1_it};
    }
    CrossIt &operator++() noexcept {
      m_p1 = m_next_p1;
      m_next_p1 = *m_prev_row1_it;

      ++m_prev_row1_it;
      ++m_prev_row2_it;
      return *this;
    }
    bool operator==(const CrossIt &other) {
      return m_prev_row1_it == other.m_prev_row1_it;
    }
  };

  const FuncT m_t;
  const FuncX m_x;
  const FuncF m_f;

  GridShape m_grid_shape;
  Comm m_comm;

public:
  TransEqSolver(const FuncT &t, const FuncX &x, const FuncF &f,
                const GridShape &grid_shape, const Comm &comm)
      : m_t(t), m_x(x), m_f(f), m_grid_shape(grid_shape), m_comm(comm) {}

private:
  Float angle_method_step(Float u1, Float u2, Float f_val) const noexcept {
    return u2 + m_grid_shape.dt() * (f_val + (u1 - u2) / m_grid_shape.dx());
  }
  Float cross_method_step(Float u1, Float u2, Float u3,
                          Float f_val) const noexcept {
    return u2 + m_grid_shape.dt() * (f_val * 2 + (u1 - u3) / m_grid_shape.dx());
  }

  template <class AngleIt, class OutIt>
  void compute_row_with_angle_method(AngleIt it, AngleIt end, OutIt out_it) {
    for (; it != end; ++it, ++out_it) {
      AnglePoints p = *it;

      *out_it = angle_method_step(p[0].u, p[1].u, m_f(p[0].t, p[1].x));
    }
  }

  template <class CrossIt, class OutIt>
  void compute_row_with_cross_method(CrossIt it, CrossIt end, OutIt out_it) {
    for (; it != end; ++it, ++out_it) {
      CrossPoints p = *it;

      *out_it = cross_method_step(p[0].u, p[1].u, p[2].u, m_f(p[0].t, p[1].x));
    }
  }

public:
  std::vector<Float> solve() {
    std::vector<Float> grid_data{};
    size_t total_grid_size = m_grid_shape.t_size() * m_grid_shape.x_size();
    grid_data.reserve(total_grid_size);

    size_t x_size = m_grid_shape.x_size();
    size_t t_size = m_grid_shape.t_size();

    // Compute first row - with given border function
    for (size_t x_idx = 0; x_idx != x_size; ++x_idx) {
      Float x_coord = m_grid_shape.x_coord(x_idx);
      grid_data.push_back(m_x(x_coord));
    }

    m_comm.set_curr_last(grid_data.back());

    // Compute other rows
    for (size_t t_idx = 1; t_idx != t_size; ++t_idx) {
      Float t_coord = m_grid_shape.t_coord(t_idx);
      Float prev_t_coord = m_grid_shape.t_coord(t_idx - 1);
      // Prepare iterators
      RowIt row1_it(grid_data.begin(), m_grid_shape, t_idx - 1);
      RowIt row1_end(grid_data.begin(), m_grid_shape, t_idx);

      grid_data.push_back(std::nan(""));
      Float &row_first = grid_data.back();

      auto curr_row_it = std::inserter(grid_data, grid_data.end());

      if (t_idx == 1) {
        AngleIt it{row1_it};
        AngleIt end{row1_end};

        compute_row_with_angle_method(it, end, curr_row_it);

        Float u1 = m_comm.get_prev_last();
        if (std::isnan(u1)) {
          row_first = m_t(t_coord);
        } else {
          Float u2 = grid_data[(t_idx - 1) * x_size];
          Float f_val = m_f(prev_t_coord, m_grid_shape.x_coord(0));

          row_first = angle_method_step(u1, u2, f_val);
        }
      } else {
        RowIt row2_it(grid_data.begin(), m_grid_shape, t_idx - 2);
        RowIt row2_end = row1_it;

        CrossIt it{row1_it, row2_it};
        CrossIt end{row1_end, row2_end};

        compute_row_with_cross_method(it, end, curr_row_it);

        Float u1 = m_comm.get_prev_last();
        if (std::isnan(u1)) {
          row_first = m_t(t_coord);
        } else {
          Float u2 = grid_data[(t_idx - 2) * x_size];
          Float u3 = grid_data[(t_idx - 1) * x_size + 1];
          Float f_val = m_f(prev_t_coord, m_grid_shape.x_coord(0));

          row_first = cross_method_step(u1, u2, u3, f_val);
        }

        u1 = grid_data[t_idx * x_size - 2];
        Float f_val = m_f(prev_t_coord, m_grid_shape.x_coord(x_size - 1));
        Float u3 = m_comm.get_next_first();
        if (std::isnan(u3)) {
          Float u2 = grid_data[t_idx * x_size - 1];
          grid_data.push_back(angle_method_step(u1, u2, f_val));
        } else {
          Float u2 = grid_data[(t_idx - 1) * x_size - 1];
          grid_data.push_back(cross_method_step(u1, u2, u3, f_val));
        }
      }

      m_comm.set_curr_first(row_first);
      m_comm.set_curr_last(grid_data.back());
    }

    return grid_data;
  }
};

} // namespace diff_compute

#endif // #ifndef DIFF_SOLVE_HPP_INCL
