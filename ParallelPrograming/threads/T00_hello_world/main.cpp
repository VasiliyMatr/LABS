
#include <vector>

#include <thread>
#include <iostream>
#include <sstream>

class ThreadJob {
  const size_t m_job_idx;

public:
  ThreadJob(size_t job_idx) : m_job_idx(job_idx) {}

  void operator()() {
    std::stringstream msg;
    msg << "Job with idx " << m_job_idx << " executed" << std::endl;
    std::cout << msg.str();
  }
};

int main() {
  size_t jobs_num = 0;
  std::cin >> jobs_num;

  if (!std::cin) {
    std::cout << "Invalid inputs!" << std::endl;
    std::terminate();
  }

  std::vector<std::thread> threads{};

  for (size_t job_idx = 0; job_idx != jobs_num; ++job_idx) {
    threads.emplace_back(ThreadJob(job_idx));
  }

  for (size_t job_idx = 0; job_idx != jobs_num; ++job_idx) {
    threads[job_idx].join();
  }

  return 0;
}
