
#include <vector>

#include <iostream>
#include <thread>

class ThreadJob {
  double m_sum = 0;
  const size_t m_begin;
  const size_t m_end;
  double &m_dest;

public:
  ThreadJob(size_t begin, size_t end, double &dest)
      : m_begin(begin), m_end(end), m_dest(dest) {}

  void operator()() {
    for (size_t i = m_begin; i != m_end; ++i) {
      m_sum += 1. / (i + 1);
    }

    m_dest = m_sum;
  }
};

int main() {
  size_t jobs_num = 0;
  size_t N = 0;

  std::cin >> jobs_num;
  std::cin >> N;

  if (!std::cin || jobs_num == 0) {
    std::cout << "Invalid inputs!" << std::endl;
    std::terminate();
  }

  size_t job_step = N / jobs_num;

  std::vector<double> threads_dests(jobs_num, 0.);
  std::vector<std::thread> threads{};

  for (size_t th_idx = 0; th_idx != jobs_num - 1; ++th_idx) {
    size_t begin = th_idx * job_step;
    size_t end = begin + job_step;

    threads.emplace_back(ThreadJob(begin, end, threads_dests[th_idx]));
  }
  {
    size_t last_idx = jobs_num - 1;
    size_t last_begin = last_idx * job_step;
    size_t last_end = N;

    threads.emplace_back(
        ThreadJob(last_begin, last_end, threads_dests[last_idx]));
  }

  for (auto &&th : threads) {
    th.join();
  }

  double total_sum = 0;

  for (auto th_sum : threads_dests) {
    total_sum += th_sum;
  }

  std::cout << "Harm series first " << N << " elements sum: " << total_sum
            << std::endl;

  return 0;
}