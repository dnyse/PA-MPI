#include "mpi.h"
#include <stdio.h>

#define SIZE 10

int main(void)
{
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double a[SIZE * SIZE], b[SIZE * SIZE];
    int displs[SIZE];
    int lens[SIZE];
    MPI_Datatype tril;
    MPI_Status status;

    for (int i = 0; i < SIZE * SIZE; i++) {
        a[i] = i + 1;
        b[i] = 0.0;
    }

    for (int i = 0; i < SIZE; i++) {
        displs[i] = i * SIZE;
        lens[i]   = i + 1;
    }

    MPI_Type_indexed(SIZE, lens, displs, MPI_DOUBLE, &tril);
    MPI_Type_commit(&tril);

    MPI_Sendrecv(a, 1, tril, rank, 0, b, 1, tril, rank, 0, MPI_COMM_WORLD, &status);

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (j <= i) {
                printf("%f ", b[i * SIZE + j]);
            } else
                printf(" ");
        }
        printf("\n");
    }

    MPI_Type_free(&tril);
    MPI_Finalize();
}
