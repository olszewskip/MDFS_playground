// Copyright: Maciek M.
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int my_rank;
  MPI_Status status;
  int N = atoi(argv[1]);
  double a[N];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
    a[99] = 123.;
    MPI_Send(a, 100, MPI_DOUBLE, 1, 17, MPI_COMM_WORLD);
  } else if (my_rank == 1) {
    MPI_Recv(a, 100, MPI_DOUBLE, 0, 17, MPI_COMM_WORLD, &status);
    printf("In process %d a[99]=%f\n", my_rank, a[99]);
  }
  MPI_Finalize();
}
