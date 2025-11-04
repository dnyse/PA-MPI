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

double dmvm(double *restrict y, const double *restrict a,
            const double *restrict x, int N, int iter) {
  double ts, te;
  int size, num, rest, rank, upperNeighbor, lowerNeighbor, Nlocal, cs, Ncurrent;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  num = N / size;
  rest = N % size;
  upperNeighbor = (rank - 1) % size;
  if (upperNeighbor < 0)
    upperNeighbor = size - 1;
  lowerNeighbor = (rank + 1) % size;
  Nlocal = (N / size) + ((N % size > rank) ? 1 : 0);
  cs = rank * (N / size) + MIN(rest, rank);
  Ncurrent = Nlocal;

  for (int rot = 0; rot < size; rot++) {
    for (int r = 0; r < Nlocal; r++) {
      for (int c = cs; c < cs + Ncurrent; c++) {
        y[r] += a[r * N + c] * x[c - cs];
      }
    }
    cs += Ncurrent;
    if (cs >= N)
      cs = 0; // wrap around
    Ncurrent = (N / size) + ((N % size > (rank + rot + 1) % size) ? 1 : 0);
    if (rot != size - 1)
      MPI_Sendrecv_replace(x, num + (rest ? 1 : 0), MPI_DOUBLE, upperNeighbor,
                           0, lowerNeighbor, 0, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
  }

  ts = getTimeStamp();
#ifdef SERIAL
  for (int j = 0; j < iter; j++) {
    for (int r = 0; r < N; r++) {
      for (int c = 0; c < N; c++) {
        y[r] = y[r] + a[r * N + c] * x[c];
      }
    }
  }
#endif
#ifdef CHECK
  {
    double sum = 0.0;

    for (int i = 0; i < N; i++) {
      sum += y[i];
      y[i] = 0.0;
    }
    fprintf(stderr, "Sum: %f\n", sum);
  }
#endif
  te = getTimeStamp();
  return te - ts;
}
