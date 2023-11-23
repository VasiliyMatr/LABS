#ifndef COMMON_HPP
#define COMMON_HPP

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <omp.h>

static constexpr size_t ISIZE = 0x2000;
static constexpr size_t JSIZE = 0x2000;

static constexpr size_t IDIST = 2;
static constexpr size_t JDIST = 3;

using A = double[ISIZE][JSIZE];

A a;

void init_a(A &a) {
    for (size_t i = 0; i != ISIZE; ++i) {
        for (size_t j = 0; j != JSIZE; ++j) {
            a[i][j] = 10 * i + j;
        }
    }
}

static void dump_a([[maybe_unused]] A &a,
                   [[maybe_unused]] const char *dump_file_name) {
#if 0
    std::ofstream out{dump_file_name};
    out << std::setprecision(3);
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

    return std::sin(0.1 * x);
}

#endif // COMMON_HPP
