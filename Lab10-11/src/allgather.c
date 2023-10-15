#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
    int size, rank;    /* total numer of procs and current proc numer */
    int parcel = 0;    /* the value to be send */
    double start, end; /* start and stop the transfer time*/
    const int parent_pid = 0;
    const int tag = 31415926;
    const int recvcount = 1;

    MPI_Init(&argc, &argv);               /* starts MPI */
    MPI_Comm_size(MPI_COMM_WORLD, &size); /* get total num of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* get number of processes */

    /* generate parcel value */
    srand(time(NULL) + rank);
    parcel = rand() % 1025;

    /* We're saving all the packages we've received */
    int parcel_collection[size];
    memset(&parcel_collection, 0, size);

    printf("Process %d going to send %4d to everyone\n", rank, parcel);

    start = MPI_Wtime(); // fixing the "parcel start" time,
                         // locally for each process

    MPI_Allgather(&parcel, 1, MPI_INT, parcel_collection, recvcount, MPI_INT, MPI_COMM_WORLD);

    end = MPI_Wtime(); // fixing the "end of reception" time,
                       // locally for each process

    /* running time of the current process */
    double my_runtime = end-start;

    printf("Process %d collected the following values:", rank);
    for (int i = 0; i < size; i++)
        printf(" %4d", parcel_collection[i]);
    printf(" (runtime was %lfsec).\n", my_runtime);

    /* submit work time for processing */
    MPI_Send(&my_runtime, 1, MPI_DOUBLE, parent_pid, tag, MPI_COMM_WORLD);

    /* Only the parent process will process the total time */
    if (rank == parent_pid)
    {
        /* Revive the work time from everyone and calc the avg */
        double sum = 0;
        double tmp = 0;
        for (int i = 0; i != size; i++)
        {
            MPI_Recv(&tmp, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += tmp;
        }
        double runtime = sum / size;
    
        printf("---\nAvg runtime was %lf.\n", runtime);
    }

    /* completes MPI */
    MPI_Finalize();
    return 0;
}