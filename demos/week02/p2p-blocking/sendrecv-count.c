#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

/*  Experiments sendind multiple elements with the sendrecv communication pattern.
 *  Optional program arguments <sendCount> <recvCount>
 *  Things to try:
 *  - Run without program arguments
 *  - Run with sending less elements than space in the recv buffer (e.g., 2 4)
 *  - Run with sending more elements than space in the recv buffer (e.g., 4 2)
 *  - Run with exchanging 100, 1000, 10000, 100000, 1000000 elements. (100 100, 1000 1000,
 * ...) What happens?
 *  - Try again with another MPI implementation (module unload intel intelmpi,
 * module load gcc/14.2.0 openmpi/4.1.3-gcc14.2.0). What happens now?
 */

int main(int argc, char** argv)
{
    int rank, size;
    size_t sendCount = 1, recvCount = 1;
    int tag = 0;
    MPI_Status status;
    int *recvbuf = NULL, *sendbuf = NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size != 2) {
        printf("This application is meant to be run with 2 processes.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (argc > 2) {
        sendCount = atoi(argv[1]);
        recvCount = atoi(argv[2]);

        if (!rank) {
            printf("Allocating %lu elements in send buffer.\n", sendCount);
            printf("Allocating %lu elements in recv buffer.\n", recvCount);
        }
    }

    sendbuf = (int*)malloc(sendCount * sizeof(int));
    recvbuf = (int*)malloc(recvCount * sizeof(int));

    for (int i = 0; i < sendCount; i++) {
        sendbuf[i] = rank;
    }
    printf("I am %d of %d and send %d\n", rank, size, sendbuf[0]);

    for (int i = 0; i < recvCount; i++) {
        recvbuf[i] = rank;
    }

    if (rank == 0) {
        MPI_Send((void*)sendbuf, sendCount, MPI_INTEGER, 1, tag, MPI_COMM_WORLD);
        MPI_Recv((void*)recvbuf, recvCount, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, &status);
    } else {
        MPI_Send((void*)sendbuf, sendCount, MPI_INTEGER, 0, tag, MPI_COMM_WORLD);
        MPI_Recv((void*)recvbuf, recvCount, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &status);
    }

    int received;
    MPI_Get_count(&status, MPI_INTEGER, &received);
    printf("I am %d and received %d elements\n", rank, received);

    free(sendbuf);
    free(recvbuf);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
