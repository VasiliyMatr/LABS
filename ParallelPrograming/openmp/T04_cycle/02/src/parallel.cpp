#include <iterator>
#include <omp.h>
#include <vector>

#include "common.hpp"

A a_copy;

int main() {
    init_a(a_copy);

    auto begin_time = omp_get_wtime();

    omp_set_num_threads(6);

#pragma omp parallel for schedule(guided)
    for (size_t i = 0; i != ISIZE - IDIST; ++i) {
        for (size_t j = JDIST; j != JSIZE; ++j) {
            a[i][j] = a_f(a_copy[i + IDIST][j - JDIST]);
        }
    }

    for (size_t i = 0; i != ISIZE; ++i) {
        for (size_t j = 0; j != JDIST; ++j) {
            a[i][j] = 10 * i + j;
        }
    }

    for (size_t i = ISIZE - IDIST; i != ISIZE; ++i) {
        for (size_t j = 0; j != JSIZE; ++j) {
            a[i][j] = 10 * i + j;
        }
    }

    auto end_time = omp_get_wtime();

    std::cout << end_time - begin_time << std::endl;

    dump_a(a, "par.txt");

    return 0;
}
