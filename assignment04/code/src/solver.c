/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "allocate.h"
#include "mpi.h"
#include "parameter.h"
#include "solver.h"

#define PI 3.14159265358979323846
#define P(i, j) p[(j) * (imax + 2) + (i)]
#define RHS(i, j) rhs[(j) * (imax + 2) + (i)]
#define DEBUG

static int sizeOfRank(int rank, int size, int N) {
  // If the rank is smaller than the rest, receive an additional row
  return N / size + ((N % size > rank) ? 1 : 0);
}

static void exchange(Solver *solver) {
  MPI_Request requests[4] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                             MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  /* exchange ghost cells with top neighbor
   * (0,0) bottom left,
   * (imax + 2, jmax + 2) top right */
  if (solver->rank + 1 < solver->size) {
    int top = solver->rank + 1;
    double *src = solver->p + (solver->jmaxLocal) * (solver->imax + 2) + 1;
    double *dst = solver->p + (solver->jmaxLocal + 1) * (solver->imax + 2) + 1;

    MPI_Isend(src, solver->imax, MPI_DOUBLE, top, 0, MPI_COMM_WORLD,
              &requests[0]);
    MPI_Irecv(dst, solver->imax, MPI_DOUBLE, top, 1, MPI_COMM_WORLD,
              &requests[1]);
  }

  /* exchange ghost cells with bottom neighbor */
  if (solver->rank > 0) {
    int bottom = solver->rank - 1;
    double *src = solver->p + (solver->imax + 2) + 1;
    double *dst = solver->p + 1;

    MPI_Isend(src, solver->imax, MPI_DOUBLE, bottom, 1, MPI_COMM_WORLD,
              &requests[2]);
    MPI_Irecv(dst, solver->imax, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD,
              &requests[3]);
  }

  MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
}

void getResult(Solver *solver) {
  double *p = NULL;
  int *rcvCounts, *displs;

  if (solver->rank == 0) {
    p = allocate(64, (solver->imax + 2) * (solver->jmax + 2) * sizeof(double));
    rcvCounts = (int *)malloc((solver->size) * sizeof(int));
    displs = (int *)malloc(solver->size * sizeof(int));
    rcvCounts[0] = solver->imax * solver->jmaxLocal;
    displs[0] = 0;
    int cursor = rcvCounts[0];

    for (int i = 1; i < solver->size; i++) {
      rcvCounts[i] = solver->imax * sizeOfRank(i, solver->size, solver->jmax);
      displs[i] = cursor;
      cursor += rcvCounts[i];
    }
  }
  // if (solver->rank == 1 || solver->rank == 2) {
  //   int imax = solver->imax + 2;
  //
  //   p = solver->p;
  //   for (int j = 1; j < solver->jmaxLocal + 1; j++) {
  //     printf("P of rank %d at row %d ", solver->rank, j);
  //     for (int i = 1; i < solver->imax + 1; i++) {
  //       printf("%f ", P(i, j));
  //     }
  //     printf("\n");
  //   }
  // }

  int cnt = solver->imax * solver->jmaxLocal;
  double *sendbuffer = (double *)malloc(cnt * sizeof(double));
  int idx = 0;
  for (int j = 1; j <= solver->jmaxLocal; j++) {
    for (int i = 1; i <= solver->imax; i++) {
      sendbuffer[idx++] = solver->p[j * (solver->imax + 2) + i];
    }
  }
  MPI_Gatherv(sendbuffer, cnt, MPI_DOUBLE, p, rcvCounts, displs, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
  if (solver->rank == 0) {
    writeResult(solver, p, "p.dat");
    free(p);
  }
  free(sendbuffer);
}

void initSolver(Solver *solver, Parameter *params, int problem) {
  MPI_Comm_rank(MPI_COMM_WORLD, &(solver->rank));
  MPI_Comm_size(MPI_COMM_WORLD, &(solver->size));
  solver->imax = params->imax;
  solver->jmax = params->jmax;
  solver->jmaxLocal = sizeOfRank(solver->rank, solver->size, solver->jmax);
  printf("RANK %d: imaxLocal : %d, jmaxLocal : %d\n", solver->rank,
         solver->imax, solver->jmaxLocal);

  solver->dx = params->xlength / params->imax;
  solver->dy = params->ylength / params->jmax;
  solver->ys = solver->rank * solver->jmaxLocal * solver->dy;
  solver->eps = params->eps;
  solver->omega = params->omg;
  solver->itermax = params->itermax;

  int imax = solver->imax;
  int jmax = solver->jmax;
  int jmaxLocal = solver->jmaxLocal;

  // adapt for MPI case
  size_t bytesize = (imax + 2) * (jmaxLocal + 2) * sizeof(double);
  solver->p = allocate(64, bytesize);
  solver->rhs = allocate(64, bytesize);

  double dx = solver->dx;
  double dy = solver->dy;
  double *p = solver->p;
  double *rhs = solver->rhs;

  // adapt for MPI case
  double y_global;
  for (int j = 0; j < jmaxLocal + 2; j++) {
    for (int i = 0; i < imax + 2; i++) {
      y_global = solver->ys + j * dy;
      P(i, j) = sin(2.0 * PI * i * dx * 2.0) + sin(2.0 * PI * y_global * 2.0);
    }
  }

  if (problem == 2) {
    for (int j = 0; j < jmaxLocal + 2; j++) {
      for (int i = 0; i < imax + 2; i++) {
        RHS(i, j) = sin(2.0 * PI * i * dx);
      }
    }
  } else {
    for (int j = 0; j < jmaxLocal + 2; j++) {
      for (int i = 0; i < imax + 2; i++) {
        RHS(i, j) = 0.0;
      }
    }
  }
}

void solve(Solver *solver) {
  int imax = solver->imax;
  int jmax = solver->jmax;
  int jmaxLocal = solver->jmaxLocal;
  double eps = solver->eps;
  int itermax = solver->itermax;
  double dx2 = solver->dx * solver->dx;
  double dy2 = solver->dy * solver->dy;
  double idx2 = 1.0 / dx2;
  double idy2 = 1.0 / dy2;
  double factor = solver->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
  double *p = solver->p;
  double *rhs = solver->rhs;
  double epssq = eps * eps;
  int it = 0;
  double res = eps + 1.0;

  while ((res >= epssq) && (it < itermax)) {
    res = 0.0;
    exchange(solver);

    // adapt for mpi
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imax + 1; i++) {
        double r =
            RHS(i, j) - ((P(i - 1, j) - 2.0 * P(i, j) + P(i + 1, j)) * idx2 +
                         (P(i, j - 1) - 2.0 * P(i, j) + P(i, j + 1)) * idy2);

        P(i, j) -= (factor * r);
        res += (r * r);
      }
    }
    for (int i = 1; i < imax + 1; i++) {
      P(i, 0) = P(i, 1);
      P(i, jmaxLocal + 1) = P(i, jmaxLocal);
    }

    for (int j = 1; j < jmaxLocal + 1; j++) {
      P(0, j) = P(1, j);
      P(imax + 1, j) = P(imax, j);
    }

    MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    res = res / (double)(imax * jmax);
#ifdef DEBUG
    if (solver->rank == 0) {
      printf("%d Residuum: %e\n", it, res);
    }
#endif
    it++;
  }

  if (solver->rank == 0) {
    printf("Solver took %d iterations to reach %f using omega=%f\n", it,
           sqrt(res), solver->omega);
  }
}

void writeResult(Solver *solver, double *p_global, char *filename) {
  int imax = solver->imax;
  int jmax = solver->jmax;

  FILE *fp;
  fp = fopen(filename, "w");

  if (fp == NULL) {
    printf("Error!\n");
    exit(EXIT_FAILURE);
  }

  for (int j = 0; j < jmax; j++) {
    for (int i = 0; i < imax; i++) {
			fprintf(fp, "%f ", p_global[j * imax + i]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}
