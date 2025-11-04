#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char** argv)
{
    int count = 1;
    int tag   = 0;
    int rank, size;
    int src = 0, dst = 0;
    MPI_Status status;
    int *recvbuf, *sendbuf;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
        MPI_Finalize();
        return EXIT_SUCCESS;
    }

    sendbuf = (int*)malloc(count * sizeof(int));
    recvbuf = (int*)malloc(count * sizeof(int));

    for (int i = 0; i < count; i++) {
        sendbuf[i] = rank;
        recvbuf[i] = rank;
    }

    if (rank == 0) {
        dst = 1;
        src = 1;
    } else if (rank == 1) {
        dst = 0;
        src = 0;
    }
    printf("I am %d and send to %d: %d\n", rank, dst, sendbuf[0]);
    int received;

    for (int i = 0; i < 10; i++) {
        MPI_Sendrecv(sendbuf,
            count,
            MPI_INT,
            dst,
            tag,
            recvbuf,
            count,
            MPI_INT,
            src,
            tag,
            MPI_COMM_WORLD,
            &status);

        MPI_Get_count(&status, MPI_INT, &received);
        printf("I am %d and received %d elements\n", rank, received);
    }

    free(sendbuf);
    free(recvbuf);

    MPI_Finalize();
}
