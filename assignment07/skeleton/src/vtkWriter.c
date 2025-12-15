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
  printf("[DEBUG] Rank %d: resetFileview - syncing file\n", o->comm.rank);
  fflush(stdout);
  MPI_File_sync(o->fh);
  printf("[DEBUG] Rank %d: resetFileview - barrier\n", o->comm.rank);
  fflush(stdout);
  MPI_Barrier(o->comm.comm);
  MPI_File_get_size(o->fh, &disp);
  printf("[DEBUG] Rank %d: resetFileview - file size = %lld\n", o->comm.rank, (long long)disp);
  fflush(stdout);
  MPI_File_set_view(o->fh, disp, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  printf("[DEBUG] Rank %d: resetFileview - DONE\n", o->comm.rank);
  fflush(stdout);
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

  printf("[DEBUG] Rank %d: writeHeader - START\n", o->comm.rank);
  fflush(stdout);

  // Only rank 0 should write this header to the file.
  // Use MPI IO along with rank check to write this header to the vtk file.

  MPI_Offset disp;
  MPI_File_get_size(o->fh, &disp);
  printf("[DEBUG] Rank %d: writeHeader - file size = %lld\n", o->comm.rank, (long long)disp);
  fflush(stdout);
  
  MPI_File_set_view(o->fh, disp, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  char header[1024]; // Adjust size if needed
  int offset = 0;

  offset += snprintf(header + offset, sizeof(header) - offset,
                     "# vtk DataFile Version 3.0\n");
  offset += snprintf(header + offset, sizeof(header) - offset,
                     "PAMPI cfd solver output\n");

  printf("[DEBUG] Rank %d: writeHeader - fmt = %d\n", o->comm.rank, o->fmt);
  fflush(stdout);

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

  printf("[DEBUG] Rank %d: writeHeader - header length = %d\n", o->comm.rank, offset);
  fflush(stdout);

  if (commIsMaster(&o->comm)) {
    printf("[DEBUG] Rank %d: writeHeader - MASTER writing header\n", o->comm.rank);
    fflush(stdout);
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }
  
  printf("[DEBUG] Rank %d: writeHeader - barrier before exit\n", o->comm.rank);
  fflush(stdout);
  MPI_Barrier(o->comm.comm);
  printf("[DEBUG] Rank %d: writeHeader - DONE\n", o->comm.rank);
  fflush(stdout);
}

void vtkOpen(VtkOptions *o, char *problem) {

  printf("[DEBUG] Rank %d: vtkOpen - START for problem '%s'\n", o->comm.rank, problem);
  fflush(stdout);

  // Open the vtk file using MPI IO
  // Use MPI_File_open() in CREATE or WRONLY mode
  char filename[50];
  snprintf(filename, 50, "%s.vtk", problem);
  
  printf("[DEBUG] Rank %d: vtkOpen - opening file '%s'\n", o->comm.rank, filename);
  fflush(stdout);
  
  int err = MPI_File_open(o->comm.comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &(o->fh));
  
  if (err != MPI_SUCCESS) {
    char error_string[BUFSIZ];
    int length;
    MPI_Error_string(err, error_string, &length);
    fprintf(stderr, "[ERROR] Rank %d: MPI_File_open failed: %s\n", o->comm.rank, error_string);
    fflush(stderr);
    MPI_Abort(o->comm.comm, err);
  }
  
  printf("[DEBUG] Rank %d: vtkOpen - file opened successfully, fh=%p\n", o->comm.rank, (void*)o->fh);
  fflush(stdout);

  writeHeader(o);

  if (commIsMaster(&o->comm)) {
    printf("Writing VTK output for %s\n", problem);
  }
  
  printf("[DEBUG] Rank %d: vtkOpen - DONE\n", o->comm.rank);
  fflush(stdout);
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
  printf("[DEBUG] Rank %d: vtkScalar - START for '%s'\n", o->comm.rank, name);
  fflush(stdout);
  
  // Write MPI IO code to parallelise writing to a single file.
  // Steps to perform MPI IO
  // 1. Use resetFileview(o); here before starting MPI IO
  resetFileview(o);
  
  if (!isInitialized(&(o->fh)))
    return;

  printf("[DEBUG] Rank %d: vtkScalar - file handle initialized\n", o->comm.rank);
  fflush(stdout);

  // 2. Write header (only rank 0)
  if (commIsMaster(&(o->comm))) {
    char header[1024];
    int offset = 0;
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "SCALARS %s double 1\n", name);
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "LOOKUP_TABLE default\n");
    printf("[DEBUG] Rank %d: vtkScalar - MASTER writing scalar header\n", o->comm.rank);
    fflush(stdout);
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  printf("[DEBUG] Rank %d: vtkScalar - resetting fileview\n", o->comm.rank);
  fflush(stdout);

  // Reset file view for binary data
  resetFileview(o);

  int offsets[NDIMS];
  // 3. Get offsets from all ranks
  commGetOffsets(&(o->comm), offsets, o->p.imax, o->p.jmax, o->p.kmax);
  
  printf("[DEBUG] Rank %d: vtkScalar - offsets = [%d, %d, %d]\n", 
         o->comm.rank, offsets[0], offsets[1], offsets[2]);
  printf("[DEBUG] Rank %d: vtkScalar - local dims = [%d, %d, %d]\n",
         o->comm.rank, o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal);
  fflush(stdout);

  MPI_Offset disp;
  MPI_Datatype fileViewType;
  // 4. Get current file size for displacement
  MPI_File_get_size(o->fh, &disp);
  printf("[DEBUG] Rank %d: vtkScalar - displacement = %lld\n", o->comm.rank, (long long)disp);
  fflush(stdout);

  // 5. Create subarray for file view (global problem dimensions)
  printf("[DEBUG] Rank %d: vtkScalar - creating file view subarray\n", o->comm.rank);
  fflush(stdout);
  
  MPI_Type_create_subarray(
      NDIMS, (int[NDIMS]){o->p.kmax, o->p.jmax, o->p.imax},
      (int[NDIMS]){o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal},
      offsets, MPI_ORDER_C, MPI_DOUBLE, &fileViewType);
  MPI_Type_commit(&fileViewType);

  printf("[DEBUG] Rank %d: vtkScalar - setting file view\n", o->comm.rank);
  fflush(stdout);

  // 6. Set file view for each rank
  MPI_File_set_view(o->fh, disp, MPI_DOUBLE, fileViewType, "external32",
                    MPI_INFO_NULL);

  printf("[DEBUG] Rank %d: vtkScalar - preparing to copy bulk data\n", o->comm.rank);
  fflush(stdout);

  // CRITICAL FIX: Copy data to malloc'd buffer to avoid any potential issues
  // with MPI I/O accessing memory allocated by posix_memalign
  size_t bulk_size = (size_t)o->comm.kmaxLocal * o->comm.jmaxLocal * o->comm.imaxLocal;
  printf("[DEBUG] Rank %d: vtkScalar - allocating temp buffer (%zu doubles = %zu bytes)\n", 
         o->comm.rank, bulk_size, bulk_size * sizeof(double));
  fflush(stdout);
  
  double *tmp = malloc(bulk_size * sizeof(double));
  if (tmp == NULL) {
    fprintf(stderr, "[ERROR] Rank %d: Failed to allocate temp buffer\n", o->comm.rank);
    fflush(stderr);
    MPI_Abort(o->comm.comm, 1);
  }
  
  printf("[DEBUG] Rank %d: vtkScalar - malloc successful tmp=%p, copying bulk data from s=%p\n", 
         o->comm.rank, (void*)tmp, (void*)s);
  fflush(stdout);
  
  // Copy bulk data (excluding ghost cells) to contiguous buffer
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

  printf("[DEBUG] Rank %d: vtkScalar - copied %d values (expected %zu)\n", 
         o->comm.rank, idx, bulk_size);
  fflush(stdout);

  if (idx != (int)bulk_size) {
    fprintf(stderr, "[ERROR] Rank %d: Copy mismatch! idx=%d, expected=%zu\n",
            o->comm.rank, idx, bulk_size);
    fflush(stderr);
  }

  printf("[DEBUG] Rank %d: vtkScalar - writing %zu doubles from tmp buffer\n", 
         o->comm.rank, bulk_size);
  fflush(stdout);

  // 8. Write data (each rank writes its portion) - using contiguous malloc'd buffer
  MPI_File_write_all(o->fh, tmp, bulk_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

  printf("[DEBUG] Rank %d: vtkScalar - write complete, freeing tmp\n", o->comm.rank);
  fflush(stdout);
  
  free(tmp);

  printf("[DEBUG] Rank %d: vtkScalar - tmp freed, freeing file view type\n", o->comm.rank);
  fflush(stdout);

  // 9. Clean up datatypes
  MPI_Type_free(&fileViewType);

  printf("[DEBUG] Rank %d: vtkScalar - writing newline\n", o->comm.rank);
  fflush(stdout);

  // 10. Binary segment must be terminated with newline character
  resetFileview(o);
  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
  }
  
  printf("[DEBUG] Rank %d: vtkScalar - DONE\n", o->comm.rank);
  fflush(stdout);
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
  printf("[DEBUG] Rank %d: vtkVector - START for '%s'\n", o->comm.rank, name);
  fflush(stdout);
  
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

  printf("[DEBUG] Rank %d: vtkVector - local=[%d,%d,%d], global=[%d,%d,%d]\n",
         o->comm.rank, kmaxLocal, jmaxLocal, imaxLocal, kmax, jmax, imax);
  fflush(stdout);

  if (commIsMaster(&o->comm)) {
    printf("Register vector %s\n", name);
  }
  if (!isInitialized(&(o->fh)))
    return;

  printf("[DEBUG] Rank %d: vtkVector - file handle initialized\n", o->comm.rank);
  fflush(stdout);

  // 2. Write header (only rank 0)
  if (commIsMaster(&o->comm)) {
    char header[1024];
    int offset = 0;
    offset += snprintf(header + offset, sizeof(header) - offset,
                       "VECTORS %s double\n", name);
    printf("[DEBUG] Rank %d: vtkVector - MASTER writing vector header\n", o->comm.rank);
    fflush(stdout);
    MPI_File_write(o->fh, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  printf("[DEBUG] Rank %d: vtkVector - resetting fileview\n", o->comm.rank);
  fflush(stdout);

  // Reset file view for binary data
  resetFileview(o);

  // 3. Get offsets from all ranks
  int offsets[NDIMS];
  commGetOffsets(&(o->comm), offsets, imax, jmax, kmax);
  
  printf("[DEBUG] Rank %d: vtkVector - offsets = [%d, %d, %d]\n", 
         o->comm.rank, offsets[0], offsets[1], offsets[2]);
  fflush(stdout);

  // 4. Get current file size for displacement
  MPI_Offset disp;
  MPI_Datatype fileViewType, vectorType;
  MPI_File_get_size(o->fh, &disp);
  
  printf("[DEBUG] Rank %d: vtkVector - displacement = %lld\n", o->comm.rank, (long long)disp);
  fflush(stdout);

  // 5. Create a contiguous type for a 3D vector (3 doubles)
  printf("[DEBUG] Rank %d: vtkVector - creating vector type\n", o->comm.rank);
  fflush(stdout);
  
  MPI_Type_contiguous(NDIMS, MPI_DOUBLE, &vectorType);
  MPI_Type_commit(&vectorType);

  // 6. Create subarray for file view (each element is now a vector)
  printf("[DEBUG] Rank %d: vtkVector - creating file view subarray\n", o->comm.rank);
  fflush(stdout);
  
  MPI_Type_create_subarray(NDIMS, (int[NDIMS]){kmax, jmax, imax},
                           (int[NDIMS]){kmaxLocal, jmaxLocal, imaxLocal},
                           offsets, MPI_ORDER_C, vectorType, &fileViewType);
  MPI_Type_commit(&fileViewType);

  // 7. Set file view
  printf("[DEBUG] Rank %d: vtkVector - setting file view\n", o->comm.rank);
  fflush(stdout);
  
  MPI_File_set_view(o->fh, disp, vectorType, fileViewType, "external32",
                    MPI_INFO_NULL);

  // 8. Prepare data: interleave velocity components and average at cell centers
  size_t cnt = (size_t)imaxLocal * jmaxLocal * kmaxLocal;
  size_t bytes_needed = cnt * NDIMS * sizeof(double);
  
  printf("[DEBUG] Rank %d: vtkVector - allocating %zu bytes for tmp (cnt=%zu)\n", 
         o->comm.rank, bytes_needed, cnt);
  fflush(stdout);
  
  double *tmp = malloc(bytes_needed);
  if (tmp == NULL) {
    fprintf(stderr, "[ERROR] Rank %d: Failed to allocate %zu bytes for tmp\n", 
            o->comm.rank, bytes_needed);
    fflush(stderr);
    MPI_Abort(o->comm.comm, 1);
  }
  
  printf("[DEBUG] Rank %d: vtkVector - malloc successful, tmp=%p\n", o->comm.rank, (void*)tmp);
  fflush(stdout);
  
  int idx = 0;

  printf("[DEBUG] Rank %d: vtkVector - filling tmp array (vec.u=%p, vec.v=%p, vec.w=%p)\n",
         o->comm.rank, (void*)vec.u, (void*)vec.v, (void*)vec.w);
  fflush(stdout);

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

  printf("[DEBUG] Rank %d: vtkVector - filled %d values (expected %zu)\n", 
         o->comm.rank, idx, cnt * NDIMS);
  fflush(stdout);
  
  if (idx != (int)(cnt * NDIMS)) {
    fprintf(stderr, "[ERROR] Rank %d: Index mismatch! idx=%d, expected=%zu\n",
            o->comm.rank, idx, cnt * NDIMS);
    fflush(stderr);
  }

  // 9. Write vector data
  printf("[DEBUG] Rank %d: vtkVector - writing %zu vectors\n", o->comm.rank, cnt);
  fflush(stdout);
  
  MPI_File_write_all(o->fh, tmp, cnt, vectorType, MPI_STATUS_IGNORE);
  
  printf("[DEBUG] Rank %d: vtkVector - write complete, freeing tmp\n", o->comm.rank);
  fflush(stdout);

  // 10. Free temporary array and datatypes
  free(tmp);
  
  printf("[DEBUG] Rank %d: vtkVector - tmp freed, freeing MPI types\n", o->comm.rank);
  fflush(stdout);
  
  MPI_Type_free(&fileViewType);
  MPI_Type_free(&vectorType);

  printf("[DEBUG] Rank %d: vtkVector - writing newline\n", o->comm.rank);
  fflush(stdout);

  // 11. Binary segment must be terminated with newline character
  resetFileview(o);
  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
  }
  
  printf("[DEBUG] Rank %d: vtkVector - DONE\n", o->comm.rank);
  fflush(stdout);
}

void vtkClose(VtkOptions *o) {
  printf("[DEBUG] Rank %d: vtkClose - START\n", o->comm.rank);
  fflush(stdout);
  
  printf("[DEBUG] Rank %d: vtkClose - closing file handle\n", o->comm.rank);
  fflush(stdout);
  
  MPI_File_close(&o->fh);
  o->fh = NULL;
  
  printf("[DEBUG] Rank %d: vtkClose - DONE\n", o->comm.rank);
  fflush(stdout);
}
