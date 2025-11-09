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
// #define SERIAL
#define PARALLEL
// #define CHECK
#define BLOCKING
// #define NON_BLOCKING

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
	int Nmax = (N / size) + ((N % size > 0) ? 1 : 0);

  int x_start = rank * num + MIN(rest, rank);

  double *x_buffers[2]; // Two buffers for rotation
  x_buffers[0] = (double *)malloc(Nmax * sizeof(double));
#ifdef NON_BLOCKING
  x_buffers[1] = (double *)malloc(Nmax * sizeof(double));
  MPI_Request requests[2];
#endif
  ts = getTimeStamp();
  int b_idx = 0;
  for (int j = 0; j < iter; j++) {
    b_idx = 0;
    for (int i = 0; i < Nlocal; i++) {
      x_buffers[b_idx][i] = x[x_start + i];
    }

    upperNeighbor = (rank - 1) % size;
    lowerNeighbor = (rank + 1) % size;
    if (upperNeighbor < 0)
      upperNeighbor = size - 1;

    cs = x_start;
    Ncurrent = Nlocal;

    for (int rot = 0; rot < size; rot++) {
#ifdef NON_BLOCKING
      if (rot != size - 1) {
        MPI_Isend(x_buffers[b_idx], Nmax, MPI_DOUBLE, upperNeighbor, 0,
                  MPI_COMM_WORLD, &requests[0]);

        MPI_Irecv(x_buffers[(b_idx + 1) % 2], Nmax, MPI_DOUBLE, lowerNeighbor,
                  0, MPI_COMM_WORLD, &requests[1]);
      }
#endif
      for (int r = 0; r < Nlocal; r++) {
        int r_local = x_start + r;
        for (int c = cs; c < cs + Ncurrent; c++) {
          y[r_local] += a[r_local * N + c] * x_buffers[b_idx][c - cs];
        }
      }
      cs += Ncurrent;
      if (cs >= N)
        cs = 0;
      Ncurrent = (N / size) + ((N % size > (rank + rot + 1) % size) ? 1 : 0);
#ifdef NON_BLOCKING
      if (rot != size - 1)
        MPI_Waitall(2, requests, MPI_STATUS_IGNORE);
      b_idx = (b_idx + 1) % 2;
#endif
#ifdef BLOCKING
      if (rot != size - 1)
        MPI_Sendrecv_replace(x_buffers[b_idx], Nmax, MPI_DOUBLE, upperNeighbor,
                             0, lowerNeighbor, 0, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
#endif /* ifdef BLOCKING */
    }
  }
  te = getTimeStamp();

  free(x_buffers[0]);
#ifdef NON_BLOCKING
  free(x_buffers[1]);
#endif

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
