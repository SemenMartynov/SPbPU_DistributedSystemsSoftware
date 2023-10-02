#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char **argv)
{
    int size = 0;
    int rank = 0;
    int bid = 0;
    const int parent_pid = 0;
    const int tag = 31415926;
    MPI_Status status;

    MPI_Init(&argc, &argv);               /* starts MPI */
    MPI_Comm_size(MPI_COMM_WORLD, &size); /* get current process id */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* get number of processes */

    if (rank == parent_pid)
    { /* parent process */
        printf("R");

        int max_bid_val = 0;
        int max_bid_id = 0;
        for (int i = 1; i < size; ++i)
        {
            MPI_Recv(&bid, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            printf("Received bid %d from the auction participant %d\n", bid, i);

            /* update the winner */
            if (bid > max_bid_val)
            {
                max_bid_val = bid;
                max_bid_id = i;
            }
        }

        printf("Max bid is %d, the participant %d win!\n", max_bid_val, max_bid_id);
    }
    else
    { /* child process */
        srand(time(NULL) + rank);
        bid = rand() % 1025;
        printf("The auction participant %d bids %d\n", rank, bid);
        MPI_Send(&bid, 1, MPI_INT, parent_pid, tag, MPI_COMM_WORLD);
    }

    /* completes MPI */
    MPI_Finalize();
    return EXIT_SUCCESS;
}
