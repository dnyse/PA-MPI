#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char** argv)
{
    int count = 1;
    int tag   = 0;
    int rank, size;
    MPI_Status status;
    int *recvbuf, *sendbuf;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size != 2) {
        printf("This application is meant to be run with 2 processes.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (argc > 1) {
        count = atoi(argv[1]);
        if (!rank) {
            printf("Allocate %d elements.\n", count);
        }
    } else {
        if (!rank) {
            printf("Program expects buffer size as argument! "
                   "Aborting!\n");
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    sendbuf = (int*)malloc(count * sizeof(int));
    recvbuf = (int*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        sendbuf[i] = rank;
        recvbuf[i] = rank;
    }
    printf("I am %d and send %d\n", rank, sendbuf[0]);

    int peer = 1 - rank;
    int received;

    for (int i = 0; i < 10; i++) {
        MPI_Send((void*)sendbuf, count, MPI_INT, peer, tag, MPI_COMM_WORLD);
        MPI_Recv((void*)recvbuf, count, MPI_INT, peer, tag, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &received);
        printf("%d - I am %d and received %d elements\n", i, rank, received);
    }

    free(sendbuf);
    free(recvbuf);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
