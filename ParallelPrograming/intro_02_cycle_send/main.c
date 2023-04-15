#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

int main(int argc, char **argv) {
  int commsize = 0;
  int my_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int counter = 1;

  if (my_rank == 0) {
    printf("Rank 0: Send counter = %d\n", counter);
    MPI_Send((void *)&counter, 1, MPI_INT, (commsize == 1) ? 0 : 1, 0, MPI_COMM_WORLD);
    MPI_Recv(&counter, 1, MPI_INT, commsize - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Rank 0: Recv counter = %d\n", counter);
  }
  else {
    MPI_Recv(&counter, 1, MPI_INT, my_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Rank %d: Resend counter = %d + 1\n", my_rank, counter);
    ++counter;
    int next_rank = my_rank + 1;
    MPI_Send((void *)&counter, 1, MPI_INT, next_rank % commsize, 0, MPI_COMM_WORLD); 
  }

  MPI_Finalize();

  return 0;
}

