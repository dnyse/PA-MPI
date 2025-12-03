#include <mpi.h>
#include <sched.h>
#include <stdio.h>
#include <unistd.h>

#define MAX_NUM_THREADS 256

static int getProcessorID(cpu_set_t *cpu_set) {
  int processorId;

  for (processorId = 0; processorId < MAX_NUM_THREADS; processorId++) {
    if (CPU_ISSET(processorId, cpu_set)) {
      break;
    }
  }
  return processorId;
}

int affinity_getProcessorId(void) {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  sched_getaffinity(getpid(), sizeof(cpu_set_t), &cpu_set);

  return getProcessorID(&cpu_set);
}

int main(void) {
  int rank = 0, size = 0;
  char hostname[1024];
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  sleep(5);
  int pinnedId = affinity_getProcessorId();
  int runningId = sched_getcpu();

  for (int i = 0; i < size; i++) {
    if (i == rank) {
      printf("Rank %d on %s cpu set %d running on %d\n", rank, hostname,
             pinnedId, runningId);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  sleep(2);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
}
