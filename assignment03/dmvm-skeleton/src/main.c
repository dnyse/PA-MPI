/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "allocate.h"
#include "timing.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// NOTE: This file assumes either BLOCKING or NON_BLOCKING is defined
// in your Makefile or build command.

extern double dmvm(double *restrict y, const double *restrict a,
                   const double *restrict x, int N, int iter, int Nlocal,
                   int x_start, int size, int rank);

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  size_t bytesPerWord = sizeof(double);
  size_t N = 0;
  size_t iter = 1;
  double *a, *x, *y;
  double t0, t1;
  double walltime;

  if (argc > 2) {
    N = atoi(argv[1]);
    iter = atoi(argv[2]);
  } else {
    printf("Usage: %s <N> <iter>\n", argv[0]);
    exit(EXIT_SUCCESS);
  }

  int num = N / size;
  int rest = N % size;
  int Nlocal = (N / size) + ((N % size > rank) ? 1 : 0);
  int x_start = rank * num + MIN(rest, rank);

  a = (double *)allocate(ARRAY_ALIGNMENT, N * Nlocal * bytesPerWord);
  x = (double *)allocate(ARRAY_ALIGNMENT, (Nlocal + 1) * bytesPerWord);
  y = (double *)allocate(ARRAY_ALIGNMENT, (Nlocal + 1) * bytesPerWord);

  for (int i = 0; i < Nlocal; i++) {
    x[i] = (double)(x_start + i);
    y[i] = 0.0;
    for (int j = 0; j < N; j++) {
      a[i * N + j] = (double)(j + x_start + i);
    }
  }

  walltime = dmvm(y, a, x, N, iter, Nlocal, x_start, size, rank);
  if (rank == 0) {
    double flops = (double)2.0 * N * N * iter;
    printf("%zu %zu %.2f %.2f\n", iter, N, 1.0E-06 * flops / walltime,
           walltime);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
