#include <iostream>

#include <omp.h>
#include <sstream>

int main() {
#pragma omp parallel
    {
        std::ostringstream out;

        out << "Hello from thread #" << omp_get_thread_num() << "!"
                  << std::endl;

        std::cout << out.str();
    }

    return 0;
}
