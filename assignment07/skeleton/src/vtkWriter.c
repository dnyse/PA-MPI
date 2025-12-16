/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "allocate.h"
#include "comm.h"
#include "vtkWriter.h"

#define G(v, i, j, k) \
  v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

static void resetFileview(VtkOptions *o) {
  MPI_Offset disp;
  MPI_File_sync(o->fh);
  MPI_Barrier(o->comm.comm);
  MPI_File_get_size(o->fh, &disp);
  MPI_File_set_view(o->fh, disp, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
}

static void writeHeader(VtkOptions *o) {
  MPI_File_set_view(o->fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  
  char header[1024];
  int offset = 0;

  offset += snprintf(header + offset, sizeof(header) - offset,
                     "# vtk DataFile Version 3.0\n");
  offset += snprintf(header + offset, sizeof(header) - offset,
                     "PAMPI cfd solver output\n");
  offset += snprintf(header + offset, sizeof(header) - offset, "BINARY\n");
  offset += snprintf(header + offset, sizeof(header) - offset,
                     "DATASET STRUCTURED_POINTS\n");
  offset += snprintf(header + offset, sizeof(header) - offset,
                     "DIMENSIONS %d %d %d\n", o->grid.imax, o->grid.jmax,
                     o->grid.kmax);
  offset += snprintf(header + offset, sizeof(header) - offset, 
                     "ORIGIN %f %f %f\n",
                     o->grid.dx * 0.5, o->grid.dy * 0.5, o->grid.dz * 0.5);
  offset += snprintf(header + offset, sizeof(header) - offset,
                     "SPACING %f %f %f\n", o->grid.dx, o->grid.dy, o->grid.dz);
  offset += snprintf(header + offset, sizeof(header) - offset, 
                     "POINT_DATA %d\n",
                     o->grid.imax * o->grid.jmax * o->grid.kmax);

  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }
  
  MPI_Barrier(o->comm.comm);
}

void vtkOpen(VtkOptions *o, char *problem) {
  char filename[50];
  snprintf(filename, 50, "%s.vtk", problem);
  
  int err = MPI_File_open(o->comm.comm, filename, 
                          MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &(o->fh));
  
  if (err != MPI_SUCCESS) {
    char error_string[BUFSIZ];
    int length;
    MPI_Error_string(err, error_string, &length);
    fprintf(stderr, "Rank %d: MPI_File_open failed: %s\n", o->comm.rank, error_string);
    MPI_Abort(o->comm.comm, err);
  }

  writeHeader(o);

  if (commIsMaster(&o->comm)) {
    printf("Writing VTK output for %s\n", problem);
  }
}

static bool isInitialized(MPI_File *ptr) {
  if (*ptr == MPI_FILE_NULL) {
    printf("vtkWriter not initialize! Call vtkOpen first!\n");
    return false;
  }
  return true;
}

void vtkScalar(VtkOptions *o, char *name, double *s) {
  resetFileview(o);
  if (!isInitialized(&(o->fh)))
    return;

  if (commIsMaster(&(o->comm))) {
    char header[1024];
    int offset = 0;
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "SCALARS %s double 1\n", name);
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "LOOKUP_TABLE default\n");
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  resetFileview(o);

  int offsets[NDIMS];
  commGetOffsets(&(o->comm), offsets, o->p.imax, o->p.jmax, o->p.kmax);

  MPI_Offset disp;
  MPI_Datatype fileViewType;
  MPI_File_get_size(o->fh, &disp);

  MPI_Type_create_subarray(
      NDIMS, (int[NDIMS]){o->p.kmax, o->p.jmax, o->p.imax},
      (int[NDIMS]){o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal},
      offsets, MPI_ORDER_C, MPI_DOUBLE, &fileViewType);
  MPI_Type_commit(&fileViewType);

  MPI_File_set_view(o->fh, disp, MPI_DOUBLE, fileViewType, "external32",
                    MPI_INFO_NULL);

  size_t bulk_size = (size_t)o->comm.kmaxLocal * o->comm.jmaxLocal * o->comm.imaxLocal;
  double *tmp = malloc(bulk_size * sizeof(double));
  if (tmp == NULL) {
    fprintf(stderr, "Rank %d: Failed to allocate temp buffer\n", o->comm.rank);
    MPI_Abort(o->comm.comm, 1);
  }
  
  int idx = 0;
  for (int k = 1; k < o->comm.kmaxLocal + 1; k++) {
    for (int j = 1; j < o->comm.jmaxLocal + 1; j++) {
      for (int i = 1; i < o->comm.imaxLocal + 1; i++) {
        int src_idx = k * (o->comm.imaxLocal + 2) * (o->comm.jmaxLocal + 2) + 
                      j * (o->comm.imaxLocal + 2) + i;
        tmp[idx++] = s[src_idx];
      }
    }
  }

  MPI_File_write_all(o->fh, tmp, bulk_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
  free(tmp);

  MPI_Type_free(&fileViewType);

  resetFileview(o);
  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
  }
}

void vtkVector(VtkOptions *o, char *name, VtkVector vec) {
  resetFileview(o);

  int imaxLocal = o->comm.imaxLocal;
  int jmaxLocal = o->comm.jmaxLocal;
  int kmaxLocal = o->comm.kmaxLocal;
  int kmax = o->p.kmax;
  int jmax = o->p.jmax;
  int imax = o->p.imax;

  if (commIsMaster(&o->comm)) {
    printf("Register vector %s\n", name);
  }
  if (!isInitialized(&(o->fh)))
    return;

  if (commIsMaster(&o->comm)) {
    char header[1024];
    int offset = 0;
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "VECTORS %s double\n", name);
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  resetFileview(o);

  int offsets[NDIMS];
  commGetOffsets(&(o->comm), offsets, imax, jmax, kmax);

  MPI_Offset disp;
  MPI_Datatype fileViewType, vectorType;
  MPI_File_get_size(o->fh, &disp);

  MPI_Type_contiguous(NDIMS, MPI_DOUBLE, &vectorType);
  MPI_Type_commit(&vectorType);

  MPI_Type_create_subarray(NDIMS, (int[NDIMS]){kmax, jmax, imax},
                           (int[NDIMS]){kmaxLocal, jmaxLocal, imaxLocal},
                           offsets, MPI_ORDER_C, vectorType, &fileViewType);
  MPI_Type_commit(&fileViewType);

  MPI_File_set_view(o->fh, disp, vectorType, fileViewType, "external32",
                    MPI_INFO_NULL);

  size_t cnt = (size_t)imaxLocal * jmaxLocal * kmaxLocal;
  double *tmp = malloc(cnt * NDIMS * sizeof(double));
  if (tmp == NULL) {
    fprintf(stderr, "Rank %d: Failed to allocate tmp array\n", o->comm.rank);
    MPI_Abort(o->comm.comm, 1);
  }
  
  int idx = 0;
  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        // Average staggered velocities to cell centers
        tmp[idx++] = (G(vec.u, i, j, k) + G(vec.u, i - 1, j, k)) / 2.0;
        tmp[idx++] = (G(vec.v, i, j, k) + G(vec.v, i, j - 1, k)) / 2.0;
        tmp[idx++] = (G(vec.w, i, j, k) + G(vec.w, i, j, k - 1)) / 2.0;
      }
    }
  }

  MPI_File_write_all(o->fh, tmp, cnt, vectorType, MPI_STATUS_IGNORE);

  free(tmp);
  MPI_Type_free(&fileViewType);
  MPI_Type_free(&vectorType);

  resetFileview(o);
  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
  }
}

void vtkClose(VtkOptions *o) {
  MPI_File_close(&o->fh);
  o->fh = NULL;
}
