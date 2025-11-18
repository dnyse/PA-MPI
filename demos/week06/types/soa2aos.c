#include <mpi.h>
#include <stdio.h>

#define SIZE 10
#define DIM  3

int main(void)
{
    MPI_Init(NULL, NULL);
    int rank, size;
    double starttime, stoptime;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    starttime = MPI_Wtime();

    if (rank == 0) {
        double vec[SIZE * DIM];
        MPI_Datatype veC;

        for (int i = 0; i < SIZE; i++) {
            vec[(i * DIM)]     = (i * 10) + 1;
            vec[(i * DIM) + 1] = (i * 10) + 2;
            vec[(i * DIM) + 2] = (i * 10) + 3;
        }

        printf("AOS:\n");
        for (int i = 0; i < SIZE * DIM; i++) {
            printf("%02.0f ", vec[i]);
        }
        printf("\n");

        MPI_Type_vector(SIZE, 1, DIM, MPI_DOUBLE, &veC);
        MPI_Type_commit(&veC);
        MPI_Send(vec, 1, veC, 1, 0, MPI_COMM_WORLD);
        MPI_Send(vec + 1, 1, veC, 1, 0, MPI_COMM_WORLD);
        MPI_Send(vec + 2, 1, veC, 1, 0, MPI_COMM_WORLD);
        MPI_Type_free(&veC);
    }

    if (rank == 1) {
        double vecA[SIZE], vecB[SIZE], vecC[SIZE];
        MPI_Status status;
        MPI_Recv(vecA, SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(vecB, SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(vecC, SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        //
        printf("SOA:\n");
        for (int i = 0; i < SIZE; i++) {
            printf("%02.0f ", vecA[i]);
        }
        printf("\n");
        for (int i = 0; i < SIZE; i++) {
            printf("%02.0f ", vecB[i]);
        }
        printf("\n");
        for (int i = 0; i < SIZE; i++) {
            printf("%02.0f ", vecC[i]);
        }
        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    stoptime = MPI_Wtime();
    if (rank == 0) {
        printf("Took %02f s\n", stoptime - starttime);
    }
    MPI_Finalize();
}
