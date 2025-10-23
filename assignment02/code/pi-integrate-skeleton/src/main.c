/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "timing.h"

#define N 10000

double integrate(double, double);
double f(double x);

int main(int argc, char **argv) {
  double wcs, wce;
  double Pi;

  wcs = getTimeStamp();
  Pi = integrate(0.0, 1.0);
  wce = getTimeStamp();

  printf("Pi=%.15lf in %.3lf s \n", Pi, wce - wcs);
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

  double delta_x = (b - a) / N;
  double sum = 0.0;
  double x;

  for (double i = 0; i < N; ++i) {
    x = a + (i + 0.5) * delta_x;
    sum += f(x);
  }

  return sum * delta_x * 4;
}
