
#include <cassert>
#include <cmath>

#ifndef INTEGRATE_HPP_INCL
#define INTEGRATE_HPP_INCL

namespace integrate {

using Float = double;

class FloatRange {
  Float m_min = std::nan("");
  Float m_max = std::nan("");
  Float m_len = std::nan("");

public:
  constexpr FloatRange(Float min, Float max) noexcept
      : m_min(min), m_max(max), m_len(max - min) {
    assert(min <= max && "Invalid FloatRange");
  }

  constexpr Float min() const noexcept { return m_min; }
  constexpr Float max() const noexcept { return m_max; }
  constexpr Float len() const noexcept { return m_len; }
};

template <class Func> Float integrate(const Func &func, const FloatRange &range, Float step) {
  size_t steps_num = range.len() / step;
  Float range_min = range.min();
  Float integral = 0;

  for (size_t i = 0; i != steps_num; ++i) {
    Float f1 = func(i * step + range_min);
    Float f2 = func((i + 1) * step + range_min);
    Float f_avg = (f1 + f2) / 2;

    integral += f_avg * step;
  }

  return integral;
}

} // namespace integrate

#endif // #ifndef INTEGRATE_HPP_INCL
