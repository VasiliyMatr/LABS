
#include <algorithm>
#include <numeric>
#include <vector>

#include <iostream>
#include <thread>

#include <integrate.hpp>

using integrate::Float;
using integrate::FloatRange;

static constexpr FloatRange INTEGRATE_RANGE{0.01, 10};

struct Func {
  constexpr Float operator()(Float x) const noexcept { return std::cos(1 / x); }
};

template <class Func> class ThreadJob {
  static constexpr Float M_BEGIN_STEP = 1E-3;
  static constexpr Float M_STEP_DIVIDER = 10;
  static constexpr Float M_MIN_STEP = 1E-10;

  const Func m_func;
  const FloatRange m_range;
  const Float m_error;
  Float &m_result;

public:
  ThreadJob(Func func, FloatRange range, Float error, Float &result)
      : m_func(func), m_range(range), m_error(error), m_result(result) {}
  void operator()() const noexcept {
    Float integral = std::nan("");
    Float prev_integral = std::nan("");

    Float step = M_BEGIN_STEP;

    for (; !(std::abs(prev_integral - integral) < m_error);
         step /= M_STEP_DIVIDER) {
      if (step < M_MIN_STEP) {
        std::cout << "WARNING: Minimal compute step reached. Can not guarantee "
                     "requested precision"
                  << std::endl;

        break;
      }

      prev_integral = integral;
      integral = integrate::integrate(m_func, m_range, step);
    }

    std::cout << "End range [" << m_range.min() << ", " << m_range.max()
              << "] with step " << step << std::endl;

    m_result = integral;
  }
};

int main() {
  size_t jobs_num = 0;
  Float error = 0;

  std::cin >> jobs_num;
  std::cin >> error;

  if (!std::cin || jobs_num == 0 || error <= 0) {
    std::cout << "Invalid inputs!" << std::endl;
    std::terminate();
  }

  Float th_error = error / jobs_num;
  Float th_range_len = INTEGRATE_RANGE.len() / jobs_num;
  std::vector<std::thread> threads;
  std::vector<Float> results(jobs_num, 0);

  for (size_t th_idx = 0; th_idx != jobs_num; ++th_idx) {
    Float th_range_min = INTEGRATE_RANGE.min() + th_range_len * th_idx;
    FloatRange th_range{th_range_min, th_range_min + th_range_len};
    threads.emplace_back(
        ThreadJob(Func(), th_range, th_error, results[th_idx]));
  }

  for (auto &&th : threads) {
    th.join();
  }

  Float integral = std::accumulate(results.cbegin(), results.cend(), Float{});

  Float digits_num = 1 - std::log10(error);
  constexpr Float DEFAULT_DIGITS_NUM = 4.;

  std::cout.precision(std::max(digits_num, DEFAULT_DIGITS_NUM));
  std::cout << "Integral = " << integral << std::endl;

  return 0;
}