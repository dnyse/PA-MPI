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

  // A: Nlocal rows (local) x N columns (global)
  a = (double *)allocate(ARRAY_ALIGNMENT, (size_t)N * Nlocal * bytesPerWord);
  // X: Full global vector (needed by everyone for initialization)
  x = (double *)allocate(ARRAY_ALIGNMENT, (size_t)N * bytesPerWord);
  // Y: Nlocal rows (local)
  y = (double *)allocate(ARRAY_ALIGNMENT, (size_t)Nlocal * bytesPerWord);

  for (int i = 0; i < N; i++) {
    x[i] = (double)i;
  }
  for (int i = 0; i < Nlocal; i++) {
    y[i] = 0.0; 
    for (int j = 0; j < N; j++) {
      a[i * N + j] = (double)j + (x_start + i);
    }
  }

  walltime = dmvm(y, a, x, N, iter, Nlocal, x_start, size, rank);
  double local_flops = (double)2.0 * N * Nlocal * iter;
  double global_flops, max_walltime;

  MPI_Reduce(&local_flops, &global_flops, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&walltime, &max_walltime, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    printf("%zu %zu %.2f %.2f\n", iter, N,
           1.0E-06 * global_flops / max_walltime, max_walltime);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
