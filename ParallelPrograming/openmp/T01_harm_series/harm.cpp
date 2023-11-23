#include <charconv>
#include <exception>
#include <iostream>
#include <numeric>
#include <system_error>
#include <vector>

#include <omp.h>

void perrAndExit(bool cond, const char *err_msg) {
    if (!cond) {
        std::cerr << err_msg;
        std::terminate();
    }
}

int main(int argc, char **argv) {
    perrAndExit(argc == 3, "Invalid arguments\n");
    size_t N = 0;
    {
        std::string n_str{argv[1]};
        auto ec =
            std::from_chars(n_str.c_str(), n_str.c_str() + n_str.size(), N).ec;

        perrAndExit(ec == std::errc(), "Invalid arguments");
        perrAndExit(N > 0, "Invalid arguments\n");
    }

    size_t num_th = 0;
    {
        std::string num_th_str{argv[2]};
        auto ec =
            std::from_chars(num_th_str.c_str(),
                            num_th_str.c_str() + num_th_str.size(), num_th)
                .ec;

        perrAndExit(ec == std::errc(), "Invalid arguments");
        perrAndExit(num_th > 0, "Invalid arguments\n");
    }

    omp_set_num_threads(num_th);

    // Avoid cache conflicts
    size_t th_padding = 8;
    std::vector<double> res_th(num_th * th_padding, 0.);

    auto start_time = omp_get_wtime();
#pragma omp parallel for
    for (size_t n = N; n != 0; --n) {
        res_th[omp_get_thread_num() * th_padding] += 1. / n;
    }

    auto res = std::accumulate(res_th.cbegin(), res_th.cend(), 0.);
    auto end_time = omp_get_wtime();

    std::cout << "H(" << N << ") = " << res << std::endl;
    std::cout << "Time = " << end_time - start_time << std::endl;

    return 0;
}
