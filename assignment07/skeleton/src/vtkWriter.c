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

#include "allocate.h"
#include "comm.h"
#include "vtkWriter.h"
// #define G(v, i, j, k) v[(k) * imax * jmax + (j) * imax + (i)]
#define G(v, i, j, k)                                                          \
  v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

// reset fileview for output of string headers
static void resetFileview(VtkOptions *o) {
  MPI_Offset disp;
  MPI_File_sync(o->fh);
  MPI_Barrier(o->comm.comm);
  MPI_File_get_size(o->fh, &disp);
  MPI_File_set_view(o->fh, disp, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
}

static double floatSwap(double f) {
  union {
    double f;
    char b[8];
  } dat1, dat2;

  dat1.f = f;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.f;
}

static void writeHeader(VtkOptions *o) {

  // Only rank 0 should write this header to the file.
  // Use MPI IO along with rank check to write this header to the vtk file.

  MPI_Offset disp;
  MPI_File_get_size(o->fh, &disp);
  MPI_File_set_view(o->fh, disp, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  char header[1024]; // Adjust size if needed
  int offset = 0;

  offset += snprintf(header + offset, sizeof(header) - offset,
                     "# vtk DataFile Version 3.0\n");
  offset += snprintf(header + offset, sizeof(header) - offset,
                     "PAMPI cfd solver output\n");

  if (o->fmt == ASCII) {
    offset += snprintf(header + offset, sizeof(header) - offset, "ASCII\n");
  } else if (o->fmt == BINARY) {
    offset += snprintf(header + offset, sizeof(header) - offset, "BINARY\n");
  }

  offset += snprintf(header + offset, sizeof(header) - offset,
                     "DATASET STRUCTURED_POINTS\n");
  offset += snprintf(header + offset, sizeof(header) - offset,
                     "DIMENSIONS %d %d %d\n", o->grid.imax, o->grid.jmax,
                     o->grid.kmax);
  offset +=
      snprintf(header + offset, sizeof(header) - offset, "ORIGIN %f %f %f\n",
               o->grid.dx * 0.5, o->grid.dy * 0.5, o->grid.dz * 0.5);
  offset += snprintf(header + offset, sizeof(header) - offset,
                     "SPACING %f %f %f\n", o->grid.dx, o->grid.dy, o->grid.dz);
  offset +=
      snprintf(header + offset, sizeof(header) - offset, "POINT_DATA %d\n",
               o->grid.imax * o->grid.jmax * o->grid.kmax);

  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(o->comm.comm);
}

void vtkOpen(VtkOptions *o, char *problem) {

  // Open the vtk file using MPI IO
  // Use MPI_File_open() in CREATE or WRONLY mode
  char filename[50];
  snprintf(filename, 50, "%s.vtk", problem);
  MPI_File_open(o->comm.comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &(o->fh));

  writeHeader(o);

  printf("Writing VTK output for %s\n", problem);
}

// static void writeScalar(VtkOptions *o, double *s) {
//
//   int imax = o->grid.imax;
//   int jmax = o->grid.jmax;
//   int kmax = o->grid.kmax;
//
//   for (int k = 0; k < kmax; k++) {
//     for (int j = 0; j < jmax; j++) {
//       for (int i = 0; i < imax; i++) {
//         if (o->fmt == ASCII) {
//           fprintf(o->fh, "%f\n", G(s, i, j, k));
//         } else if (o->fmt == BINARY) {
//           fprintf((double[1]){floatSwap(G(s, i, j, k))}, sizeof(double), 1,
//                   o->fh);
//         }
//       }
//     }
//   }
//   if (o->fmt == BINARY)
//     fprintf(o->fh, "\n");
// }
//
static bool isInitialized(MPI_File *ptr) {
  if (*ptr == MPI_FILE_NULL) {
    printf("vtkWriter not initialize! Call vtkOpen first!\n");
    return false;
  }
  return true;
}

void vtkScalar(VtkOptions *o, char *name, double *s) {
  // Write MPI IO code to parallelise writing to a single file.
  // Steps to perform MPI IO
  // 1. Use resetFileview(o); here before starting MPI IO
  resetFileview(o);
  if (!isInitialized(&(o->fh)))
    return;

  // 2. Write header (only rank 0)
  if (commIsMaster(&(o->comm))) {
    char header[1024];
    int offset = 0;
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "SCALARS %s double 1\n", name);
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "LOOKUP_TABLE default\n");
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  // Reset file view for binary data
  resetFileview(o);

  int offsets[NDIMS];
  // 3. Get offsets from all ranks
  commGetOffsets(&(o->comm), offsets, o->p.imax, o->p.jmax, o->p.kmax);

  MPI_Offset disp;
  MPI_Datatype fileViewType;
  // 4. Get current file size for displacement
  MPI_File_get_size(o->fh, &disp);

  // 5. Create subarray for file view (global problem dimensions)
  MPI_Type_create_subarray(
      NDIMS, (int[NDIMS]){o->p.kmax, o->p.jmax, o->p.imax},
      (int[NDIMS]){o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal},
      offsets, MPI_ORDER_C, MPI_DOUBLE, &fileViewType);
  MPI_Type_commit(&fileViewType);

  // 6. Set file view for each rank
  MPI_File_set_view(o->fh, disp, MPI_DOUBLE, fileViewType, "external32",
                    MPI_INFO_NULL);

  // 7. Create subarray excluding ghost cells (local data with ghosts:
  // imaxLocal+2)
  MPI_Datatype bulkType;
  MPI_Type_create_subarray(
      NDIMS,
      (int[NDIMS]){o->comm.kmaxLocal + 2, o->comm.jmaxLocal + 2,
                   o->comm.imaxLocal + 2},
      (int[NDIMS]){o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal},
      (int[NDIMS]){1, 1, 1}, MPI_ORDER_C, MPI_DOUBLE, &bulkType);
  MPI_Type_commit(&bulkType);

  // 8. Write data (each rank writes its portion)
  MPI_File_write(o->fh, s, 1, bulkType, MPI_STATUS_IGNORE);

  // 9. Clean up datatypes
  MPI_Type_free(&bulkType);
  MPI_Type_free(&fileViewType);

  // 10. Binary segment must be terminated with newline character
  resetFileview(o);
  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
  }
}

// static void writeVector(VtkOptions *o, VtkVector vec) {
//   int imax = o->grid.imax;
//   int jmax = o->grid.jmax;
//   int kmax = o->grid.kmax;
//
//   for (int k = 0; k < kmax; k++) {
//     for (int j = 0; j < jmax; j++) {
//       for (int i = 0; i < imax; i++) {
//         if (o->fmt == ASCII) {
//           fprintf(o->fh, "%f %f %f\n", G(vec.u, i, j, k), G(vec.v, i, j, k),
//                   G(vec.w, i, j, k));
//         } else if (o->fmt == BINARY) {
//           fwrite((double[3]){floatSwap(G(vec.u, i, j, k)),
//                              floatSwap(G(vec.v, i, j, k)),
//                              floatSwap(G(vec.w, i, j, k))},
//                  sizeof(double), 3, o->fh);
//         }
//       }
//     }
//   }
//   if (o->fmt == BINARY)
//     fprintf(o->fh, "\n");
// }

void vtkVector(VtkOptions *o, char *name, VtkVector vec) {
  // Write MPI IO code to parallelise writing to a single file.
  // Steps to perform MPI IO
  // 1. Use resetFileview(o); here before starting MPI IO
  resetFileview(o);

  int imaxLocal = o->comm.imaxLocal;
  int jmaxLocal = o->comm.jmaxLocal;
  int kmaxLocal = o->comm.kmaxLocal;
  int kmax = o->p.kmax;
  int jmax = o->p.jmax;
  int imax = o->p.imax;

  printf("Register vector %s\n", name);
  if (!isInitialized(&(o->fh)))
    return;

  // 2. Write header (only rank 0)
  if (commIsMaster(&o->comm)) {
    char header[1024];
    int offset = 0;
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "VECTORS %s double\n", name);
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  // Reset file view for binary data
  resetFileview(o);

  // 3. Get offsets from all ranks
  int offsets[NDIMS];
  commGetOffsets(&(o->comm), offsets, imax, jmax, kmax);

  // 4. Get current file size for displacement
  MPI_Offset disp;
  MPI_Datatype fileViewType, vectorType;
  MPI_File_get_size(o->fh, &disp);

  // 5. Create a contiguous type for a 3D vector (3 doubles)
  MPI_Type_contiguous(NDIMS, MPI_DOUBLE, &vectorType);
  MPI_Type_commit(&vectorType);

  // 6. Create subarray for file view (each element is now a vector)
  MPI_Type_create_subarray(NDIMS, (int[NDIMS]){kmax, jmax, imax},
                           (int[NDIMS]){kmaxLocal, jmaxLocal, imaxLocal},
                           offsets, MPI_ORDER_C, vectorType, &fileViewType);
  MPI_Type_commit(&fileViewType);

  // 7. Set file view
  MPI_File_set_view(o->fh, disp, vectorType, fileViewType, "external32",
                    MPI_INFO_NULL);

  // 8. Prepare data: interleave velocity components and average at cell centers
  size_t cnt = imaxLocal * jmaxLocal * kmaxLocal;
  double *tmp = allocate(64, cnt * NDIMS * sizeof(double));
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

  // 9. Write vector data
  MPI_File_write(o->fh, tmp, cnt, vectorType, MPI_STATUS_IGNORE);

  // 10. Free temporary array and datatypes
  free(tmp);
  MPI_Type_free(&fileViewType);
  MPI_Type_free(&vectorType);

  // 11. Binary segment must be terminated with newline character
  resetFileview(o);
  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
  }
}

void vtkClose(VtkOptions *o) {
  MPI_File_close(&o->fh);
  o->fh = NULL;
}
