#ifndef COMMON_HPP
#define COMMON_HPP

#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <thread>

static constexpr size_t ISIZE = 1000;
static constexpr size_t JSIZE = 1000;

static double a[ISIZE][JSIZE];
static double b[ISIZE][JSIZE];

static void init_ab() {
  for (size_t i = 0; i != ISIZE; ++i) {
    for (size_t j = 0; j != JSIZE; ++j) {
      a[i][j] = 10 * i + j;
      b[i][j] = 0;
    }
  }
}

static void dump_b([[maybe_unused]] const char *dump_file_name) {
#if 0
  std::ofstream out{dump_file_name};
  for (size_t i = 0; i != ISIZE; ++i) {
    for (size_t j = 0; j != JSIZE; ++j) {
      out << b[i][j] << " ";
    }
    out << std::endl;
  }
#endif
}

static double a_f(double x) {
  // std::this_thread::sleep_for(std::chrono::milliseconds(1));

  return std::sin(0.1 * x);
}

static double b_f(double x) {
  // std::this_thread::sleep_for(std::chrono::milliseconds(1));

  return x * 1.5;
}

#endif // COMMON_HPP
