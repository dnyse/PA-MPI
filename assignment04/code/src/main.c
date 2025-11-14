/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "likwid-marker.h"
#include "mpi.h"
#include "parameter.h"
#include "solver.h"
#include "timing.h"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  double startTime, endTime;
  Parameter params;
  Solver solver;
  initParameter(&params);

  if (argc < 2) {
    printf("Usage: %s <configFile>\n", argv[0]);
    exit(EXIT_SUCCESS);
  }
  readParameter(&params, argv[1]);
  initSolver(&solver, &params, 2);

  if (solver.rank == 0) {
    printParameter(&params);
  }

  startTime = getTimeStamp();
  solve(&solver);
  endTime = getTimeStamp();
  // writeResult(&solver, "p.dat");
  getResult(&solver);

  if (solver.rank == 0) {
    printf("Walltime %.2fs\n", endTime - startTime);
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
