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

#include "affinity.h"
#include "allocate.h"
#include "timing.h"

#define HOSTNAME_LENGTH 256

int main(int argc, char** argv)
{
    char hostname[HOSTNAME_LENGTH];
    const char* name = "World";
    if (argc > 1) {
        name = argv[1];
    }
    gethostname(hostname, HOSTNAME_LENGTH);
    printf("Hello %s from %s!\n", name, hostname);

    return EXIT_SUCCESS;
}
