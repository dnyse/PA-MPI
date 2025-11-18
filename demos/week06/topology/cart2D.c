#include "mpi.h"
#include <stdio.h>

#define SIZE 10
enum direction { LEFT = 0, RIGHT, BOTTOM, TOP, NDIRS };
enum dimension { JDIM = 0, IDIM, NDIMS };

int main(void)
{
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm comm;
    int neighbours[NDIRS];
    int coords[NDIMS];
    int dims[NDIMS]    = { 0, 0 };
    int periods[NDIMS] = { 0, 0 };
    MPI_Dims_create(size, NDIMS, dims);
    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periods, 0, &comm);
    MPI_Cart_shift(comm, IDIM, 1, &neighbours[LEFT], &neighbours[RIGHT]);
    MPI_Cart_shift(comm, JDIM, 1, &neighbours[BOTTOM], &neighbours[TOP]);
    MPI_Cart_get(comm, NDIMS, dims, periods, coords);

    fflush(stdout);
    if (!rank) {
        printf("Communication setup:\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < size; i++) {
        if (i == rank) {
            printf("\tRank %d of %d\n", rank, size);
            printf("\tNeighbours (bottom, top, left, right): %d, %d, %d, %d\n",
                neighbours[BOTTOM],
                neighbours[TOP],
                neighbours[LEFT],
                neighbours[RIGHT]);
            printf("\tCoordinates (j,i) %d, %d\n", coords[JDIM], coords[IDIM]);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (!rank) {
        int crank;
        for (int j = 0; j < dims[JDIM]; j++) {
            for (int i = 0; i < dims[IDIM]; i++) {
                MPI_Cart_rank(comm, (int[]) { j, i }, &crank);
                printf("%d->(%d,%d)\t", crank, j, i);
            }
            printf("\n");
        }
    }

    MPI_Finalize();
}
