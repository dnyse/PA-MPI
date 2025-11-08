/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include "timing.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// Define which implementation to use
#define PARALLEL
#define CHECK

double dmvm(double *restrict y, const double *restrict a,
            const double *restrict x, int N, int iter) {
  double ts, te;
#ifdef PARALLEL
  int size, num, rest, rank, upperNeighbor, lowerNeighbor, Nlocal, cs, Ncurrent;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  num = N / size;
  rest = N % size;
  Nlocal = (N / size) + ((N % size > rank) ? 1 : 0);

  int x_start = rank * num + MIN(rest, rank);

  double *x_local = (double *)malloc(Nlocal * sizeof(double));

  ts = getTimeStamp();
  for (int j = 0; j < iter; j++) { // Loop for 'iter'

    for (int i = 0; i < Nlocal; i++) {
      x_local[i] = x[x_start + i];
    }

    upperNeighbor = (rank - 1) % size;
    lowerNeighbor = (rank + 1) % size;
    if (upperNeighbor < 0)
      upperNeighbor = size - 1;

    cs = x_start;
    Ncurrent = Nlocal;

    for (int rot = 0; rot < size; rot++) {
      for (int r = 0; r < Nlocal; r++) {
        int r_local = x_start + r;
        for (int c = cs; c < cs + Ncurrent; c++) {
          y[r_local] += a[r_local * N + c] * x_local[c - cs];
        }
      }
      cs += Ncurrent;
      if (cs >= N)
        cs = 0;
      Ncurrent = (N / size) + ((N % size > (rank + rot + 1) % size) ? 1 : 0);
      if (rot != size - 1)
        MPI_Sendrecv_replace(x_local, Nlocal, MPI_DOUBLE, upperNeighbor, 0,
                             lowerNeighbor, 0, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
    }
  }
  te = getTimeStamp();

  free(x_local);

#ifdef CHECK
  double local_sum = 0.0;
  double global_sum = 0.0;

  for (int r = 0; r < Nlocal; r++) {
    int r_local = x_start + r;
    local_sum += y[r_local];
  }

  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    fprintf(stderr, "Sum: %f\n", global_sum);
  }
#endif // CHECK
#endif // PARALLEL

#ifdef SERIAL
  ts = getTimeStamp();
  for (int j = 0; j < iter; j++) {
    for (int r = 0; r < N; r++) {
      for (int c = 0; c < N; c++) {
        y[r] = y[r] + a[r * N + c] * x[c];
      }
    }
  }
  te = getTimeStamp();

#ifdef CHECK
  double sum = 0.0;
  for (int i = 0; i < N; i++) {
    sum += y[i];
  }
  fprintf(stderr, "Sum: %f\n", sum);
#endif // CHECK
#endif // SERIAL

  return te - ts;
}
