#include "common.hpp"

int main() {
    init_a(a);

    auto begin_time = omp_get_wtime();
    for (size_t i = 0; i != ISIZE - 2; ++i) {
        for (size_t j = 3; j != JSIZE; ++j) {
            a[i][j] = a_f(a[i + 2][j - 3]);
        }
    }
    auto end_time = omp_get_wtime();

    std::cout << end_time - begin_time << std::endl;

    dump_a(a, "ref.txt");

    return 0;
}
