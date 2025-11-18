#include <mpi.h>
#include <stdio.h>

#define SIZE 10
#define DIM  3

static void pack(double* in, double* out, int size, int offset)
{
    for (size_t i = 0; i < size; i++) {
        out[i] = in[(i * DIM) + offset];
    }
}

int main(void)
{
    MPI_Init(NULL, NULL);
    int rank, size;
    double starttime, stoptime;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    starttime = MPI_Wtime();

    if (rank == 0) {
        double vec[SIZE * DIM], buffer[SIZE];
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

        for (int i = 0; i < DIM; i++) {
            pack(vec, buffer, SIZE, i);
            MPI_Send(buffer, SIZE, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }
    }

    if (rank == 1) {
        double vecA[SIZE], vecB[SIZE], vecC[SIZE];
        MPI_Status status;
        MPI_Recv(vecA, SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(vecB, SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(vecC, SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

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
