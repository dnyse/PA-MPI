/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>

#include "timing.h"

double integrate(double, double);

int main (int argc, char** argv) {
    double wcs, wce;
    double  Pi;

    wcs = getTimeStamp();
    Pi = integrate(0.0, 1.0);
    wce = getTimeStamp();

    printf("Pi=%.15lf in %.3lf s \n", Pi,wce-wcs);
    return EXIT_SUCCESS;
}

double integrate(double a, double b) {
	
	/*
	
		Your logic to integrate between given interval a to b.
		Declare SLICES here and calculate delta x using a, b and SLICES.
		Iterate over number of slices, calculate the area and sum them.
		Return sum * delta x to get the value of PI.
		
	*/
	return 0.0;
}
