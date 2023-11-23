#include <iostream>

#include <omp.h>

#include "common.hpp"

int main() {
    init_a();

    auto start = omp_get_wtime();

    omp_set_num_threads(5);

#pragma omp parallel for
    for (size_t i = 0; i != ISIZE; ++i) {
        for (size_t j = 0; j != JSIZE; ++j) {
            a[i][j] = a_f(a[i][j]);
        }
    }

    auto end = omp_get_wtime();

    std::cout << "Total time: " << end - start << std::endl;

    dump_a("par.txt");

    return 0;
}
