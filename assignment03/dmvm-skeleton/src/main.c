/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

#include "allocate.h"
#include "timing.h"

extern double dmvm(double* restrict y,
    const double* restrict a,
    const double* restrict x,
    int N,
    int iter);

int main(int argc, char** argv)
{

    size_t bytesPerWord = sizeof(double);
    size_t N            = 0;
    size_t iter         = 1;
    double *a, *x, *y;
    double t0, t1;
    double walltime;

    if (argc > 2) {
        N    = atoi(argv[1]);
        iter = atoi(argv[2]);
    } else {
        printf("Usage: %s <N> <iter>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }

    a = (double*)allocate(ARRAY_ALIGNMENT, N * N * bytesPerWord);
    x = (double*)allocate(ARRAY_ALIGNMENT, N * bytesPerWord);
    y = (double*)allocate(ARRAY_ALIGNMENT, N * bytesPerWord);

    // initialize arrays
    for (int i = 0; i < N; i++) {
        x[i] = (double)i;
        y[i] = 0.0;

        for (int j = 0; j < N; j++) {
            a[i * N + j] = (double)j + i;
        }
    }

		MPI_Init(&argc, &argv);
    walltime = dmvm(y, a, x, N, iter);
		MPI_Finalize();

    double flops = (double)2.0 * N * N * iter;
    // # iterations, problem size, flop rate, walltime
    printf("%zu %zu %.2f %.2f\n", iter, N, 1.0E-06 * flops / walltime, walltime);

    return EXIT_SUCCESS;
}
