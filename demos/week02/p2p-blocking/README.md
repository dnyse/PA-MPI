# Experiments with blocking MPI P2P communication

Login to Fritz:

```shell
ssh fritz
```

Change into the folder with the demo source code.
Load the default Intel compiler and MPI modules:

```shell
module load intel intelmpi
```

Allocate one node in interactive job:

```shell
salloc -N 1 --time=01:00:00
```

## Example how to build and run program

Build source code with:

```shell
mpiicx -o ./test sendrecv.c
```

Run program with:

```shell
srun -n 2  ./test
```
