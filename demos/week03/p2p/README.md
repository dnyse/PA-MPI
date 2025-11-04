# MPI P2P communication

## `sendrecv-count.c`

- Allocate one node and load Intel MPI
  - Compile and execute the program with 2 processors
  - Increase the message size
- Allocate two nodes and load Intel MPI
  - Execute the program with 2 processors
  - Increase the message size
- How can you make the program **save** ?
- Decrease the count on the Recv of rank 1.
- Decrease the count on the Send of rank 0.
- Repeat all of the above with OpenMPI
