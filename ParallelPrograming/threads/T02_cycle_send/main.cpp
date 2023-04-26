
#include <vector>

#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>

class ThreadJob {
  static std::mutex m_mutex;
  static int m_value;

public:
  void operator()() {
    std::lock_guard<std::mutex> guard(m_mutex);
    std::cout << "Job uses value " << m_value++ << std::endl;
  }
};

std::mutex ThreadJob::m_mutex{};
int ThreadJob::m_value = 0;

int main() {
  size_t jobs_num = 0;
  std::cin >> jobs_num;

  if (!std::cin) {
    std::cout << "Invalid inputs!" << std::endl;
    std::terminate();
  }

  std::vector<std::thread> threads{};

  for (size_t job_idx = 0; job_idx != jobs_num; ++job_idx) {
    threads.emplace_back(ThreadJob{});
  }

  for (auto &&th : threads) {
    th.join();
  }

  return 0;
}
