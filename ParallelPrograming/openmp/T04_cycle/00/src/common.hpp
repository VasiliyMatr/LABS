#ifndef COMMON_HPP
#define COMMON_HPP

#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <thread>

static constexpr size_t ISIZE = 0x1000;
static constexpr size_t JSIZE = 0x1000;

static double a[ISIZE][JSIZE];

static void init_a() {
  for (size_t i = 0; i != ISIZE; ++i) {
    for (size_t j = 0; j != JSIZE; ++j) {
      a[i][j] = 10 * i + j;
    }
  }
}

static void dump_a([[maybe_unused]] const char *dump_file_name) {
#if 0
  std::ofstream out{dump_file_name};
  for (size_t i = 0; i != ISIZE; ++i) {
    for (size_t j = 0; j != JSIZE; ++j) {
      out << a[i][j] << " ";
    }
    out << std::endl;
  }
#endif
}

static double a_f(double x) {
  // std::this_thread::sleep_for(std::chrono::milliseconds(1));

  return std::sin(2 * x);
}

#endif // COMMON_HPP
