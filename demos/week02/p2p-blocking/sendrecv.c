#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

/*  Experiments with the sendrecv communication pattern.
 *  Things to try:
 *  - Run the code with 1, 2, and 3 processes. What happens? What is an simple fix, so
 * that it runs with any process count?
 *  - Change the order of the Send and Recv calls. Does the order matter?
 *  - Add more Send calls on sender and receiver side. What happens?
 *  - Add more Recv calls on sender and receiver side. What happens?
 */

int main(int argc, char** argv)
{
    int rank, size;
    int tag = 0, count = 1;
    MPI_Status status;
    int recvbuf, sendbuf;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    sendbuf = rank;
    recvbuf = rank;

    if (rank == 0) {
        MPI_Send((void*)&sendbuf, count, MPI_INTEGER, 1, tag, MPI_COMM_WORLD);
        MPI_Recv((void*)&recvbuf, count, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, &status);
    } else {
        MPI_Send((void*)&sendbuf, count, MPI_INTEGER, 0, tag, MPI_COMM_WORLD);
        MPI_Recv((void*)&recvbuf, count, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &status);
    }

    printf("I am %d of %d and received %d\n", rank, size, recvbuf);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
