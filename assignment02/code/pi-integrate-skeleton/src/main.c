/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "timing.h"

#define N 1000000000

double integrate(double, double);
double f(double x);

int main(int argc, char **argv) {
  double wcs, wce;
  double Pi, local_Pi;
  int rank, size, root_rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  double interval = 1.0 / size;

  wcs = getTimeStamp();
  local_Pi = integrate(rank * interval, (rank + 1) * interval);
  if (rank == root_rank) {
    Pi = local_Pi;

    for (int i = 1; i < size; i++) {
      double worker_result;
      MPI_Recv(&worker_result, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      Pi += worker_result;
    }
    Pi *= 4;
  } else {
    MPI_Send(&local_Pi, 1, MPI_DOUBLE, root_rank, 0, MPI_COMM_WORLD);
  }
  wce = getTimeStamp();
  if (rank == root_rank)
    printf("Pi=%.15lf in %.3lf s \n", Pi, wce - wcs);

  MPI_Finalize();
  return EXIT_SUCCESS;
}

double f(double x) { return sqrt(1 - x * x); }

double integrate(double a, double b) {

  /*

          Your logic to integrate between given interval a to b.
          Declare SLICES here and calculate delta x using a, b and SLICES.
          Iterate over number of slices, calculate the area and sum them.
          Return sum * delta x to get the value of PI.

  */
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int slices = N / size;
  double delta_x = (b - a) / slices;
  double sum = 0.0;
  double x;

  for (int i = 0; i < slices; ++i) {
    // x = a + (i + 0.5) * delta_x;
    x = a + i * delta_x;
    sum += f(x);
  }

  return sum * delta_x;
}
