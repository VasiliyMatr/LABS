#include <iostream>

#include <omp.h>
#include <sstream>
#include <thread>

int main() {
    int var = 0;

#pragma omp parallel for ordered
    for (int i = 0; i != omp_get_num_threads(); ++i) {
        std::ostringstream begin_msg;

        begin_msg << "Begin of thread #" << omp_get_thread_num() << std::endl;

        std::cout << begin_msg.str();

#pragma omp ordered
        {
            std::ostringstream var_msg;

            var += 1;
            var_msg << "Thread #" << omp_get_thread_num()
                    << " touched var = " << var << std::endl;

            std::cout << var_msg.str();
        }

        std::ostringstream end_msg;

        end_msg << "End of thread #" << omp_get_thread_num() << std::endl;

        std::cout << end_msg.str();
    }

    return 0;
}
