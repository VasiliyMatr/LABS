#include <iostream>

#include <mpi.h>

int main(int argc, char **argv) {
  int commsize = 0;
  int my_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (commsize != 2) {
    std::cout << "Invalid commsize" << std::endl;
    std::terminate();
  }

  static constexpr size_t BUFF_SIZE = 0x8;
  char buff[BUFF_SIZE] = {};
  size_t sends_num = 1000000;

  double t_begin = MPI_Wtime();

  for (size_t i = 0; i != sends_num; ++i) {
    if (my_rank == 0) {
      MPI_Send((void *)buff, sizeof(buff), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
    } else {
      MPI_Recv(buff, sizeof(buff), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }

  double t_end = MPI_Wtime();

  if (my_rank == 0) {
    std::cout << "Time per send = " << (t_end - t_begin) / sends_num
              << std::endl;
  }

  MPI_Finalize();
  return 0;
}