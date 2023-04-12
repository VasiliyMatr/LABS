#include <stdio.h>
#include <stdlib.h>

#include<mpi.h>

// Sum harmonic series elements
// begin: first element index
// end: last element index (not included in sum range)
// Example:
// sum_harm_series_elems(0,2) = 1./1 + 1./2
double sum_harm_series_elems(long long begin, long long end) {
  double sum = 0;

  for(long long i = begin; i != end; ++i) {
    sum += 1.0 / (i + 1);
  }

  return sum;
}

// Collect calculations results from all ranks except the ZERO rank
// commsize: communicator size
// result: Address to sum results to
//
// return: MPI_SUCCESS or MPI error identifier
int collect_results(int commsize, double *result) {
  for (int i = 1; i != commsize; ++i) {
    double i_result = 0;
    int error_id = MPI_Recv(&i_result, 1, MPI_DOUBLE, i, MPI_ANY_TAG,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (error_id != MPI_SUCCESS) {
      return error_id;
    }

    *result += i_result;
  }

  return MPI_SUCCESS;
}

// Send calculations result to the ZERO rank
// my_results: results to send
//
// return: MPI_SUCCESS or MPI error identifier
int send_result(double my_results) {
  return MPI_Send((void *)&my_results, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

void finalize_n_exit() {
  MPI_Finalize();
  exit(-1);
}

void perr_n_exit(int my_rank, char *msg) {
  fprintf(stderr, "Rank %d failed: %s", my_rank, msg);
  finalize_n_exit();
}

void perr_if_ZERO_n_exit(int my_rank, char *msg) {
  if (my_rank == 0) {
    perr_n_exit(my_rank, msg);
  }

  finalize_n_exit();
}

int main(int argc, char **argv) {
  int commsize = 0;
  int my_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  long long N = 0;
  if (argc != 2) {
    perr_if_ZERO_n_exit(my_rank, "Invalid arguments\n");
  }

  if (sscanf(argv[1], " %lld ", &N) != 1) {
    perr_if_ZERO_n_exit(my_rank, "Invalid arguments\n");
  }
  
  long long count_step = N / commsize;
  long long my_range_begin = count_step * my_rank;
  long long my_range_end = (my_rank == commsize - 1) ? N : my_range_begin + count_step;

  double my_range_result = sum_harm_series_elems(my_range_begin, my_range_end);

  if (my_rank == 0) {
    double result = my_range_result;

    int error_id = collect_results(commsize, &result);
    if (error_id == MPI_SUCCESS) {
      printf("Result: %lf\n", result);
    }
    else {
      perr_n_exit(my_rank, "Failed to collect results\n");
    }
  }
  else {
    int error_id = send_result(my_range_result);
    if (error_id != MPI_SUCCESS) {
      perr_n_exit(my_rank, "Failed to send results\n");
    }
  }

  MPI_Finalize();

  return 0;
}

