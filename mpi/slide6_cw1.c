// Copyright: Maciek M.
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int my_rank;
  MPI_Status status;
  int N = atoi(argv[1]);
  const int kStatistic = 100000;
  double t1, t2;
  double times[kStatistic];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int a[N];
  int i;
  for (i = 0; i < kStatistic; i++) {
    t1 = MPI_Wtime();
    if (my_rank == 0) {
      a[N - 1] = 123.;
      MPI_Send(a, N, MPI_INT, 1, 17, MPI_COMM_WORLD);
      MPI_Recv(a, N, MPI_INT, 1, 18, MPI_COMM_WORLD, &status);
    } else if (my_rank == 1) {
      MPI_Recv(a, N, MPI_INT, 0, 17, MPI_COMM_WORLD, &status);
      a[N-1] = 321;
      MPI_Send(a, N, MPI_INT, 0, 18, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();
    if (my_rank == 0) {
      double time_delta = t2 - t1;
      // printf("Send-recv in %f sec.\n", time_delta);
      times[i] = time_delta;
    }
  }
  MPI_Finalize();

  if (my_rank == 0) {
    double avg = 0, avg_sq = 0, std;
    for (i = 0; i < kStatistic; i++) {
      avg += times[i];
      avg_sq += times[i] * times[i];
    }
    avg /= kStatistic;
    avg_sq /= kStatistic;
    std = sqrt(avg_sq - avg * avg);
    printf("Avg send-recv time: %f sec. +- %f\n", avg, std);
  }
}
