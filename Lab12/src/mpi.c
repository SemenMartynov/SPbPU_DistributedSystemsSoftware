#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <float.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include "mpi.h"
#include "matrix.h"

#define PARENT_ID 0  /* taskId of first task */
#define MSG_TAG 3141 /* message tag */

int main(int argc, char **argv)
{
    int size = 0;
    int rank = 0;

    MPI_Init(&argc, &argv);               /* starts MPI */
    MPI_Comm_size(MPI_COMM_WORLD, &size); /* get current process id */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* get number of processes */

    if (rank == PARENT_ID)
    { /* parent process */

        struct Matrix m;

        /* Read input from file */
        FILE *fp = fopen("input.txt", "r");
        fscanf(fp, "%zu %zu", &m.row, &m.col);
        m.data = (double *)malloc(m.row * m.col * sizeof(double));
        for (size_t i = 0; i != m.row; ++i)
        {
            for (size_t j = 0; j != m.col; ++j)
            {
                fscanf(fp, "%lf", &m.data[i * m.col + j]);
            }
        }
        m.activeCol = m.col - 1;
        fclose(fp);

        m.passedLine = (bool *)malloc(m.row * sizeof(bool));
        memset(m.passedLine, 0, m.row * sizeof(bool));
        double main_determinant = determinant(&m, 0, m.activeCol);
        free(m.passedLine);
        if (fabs(main_determinant) < DBL_EPSILON)
        {
            printf("\nDetermining roots for a given system of equations is not possible.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        struct timespec start, finish;
        clock_gettime(CLOCK_MONOTONIC, &start);
        /* Send tasks to the workers */
        for (int pid = 1, id = 0; pid < size; pid++, id++)
        {
            MPI_Send(&m.row, 1, MPI_INT, pid, MSG_TAG, MPI_COMM_WORLD);
            MPI_Send(&m.col, 1, MPI_INT, pid, MSG_TAG, MPI_COMM_WORLD);

            MPI_Send(&(m.data[0]), m.row * m.col, MPI_DOUBLE, pid, MSG_TAG, MPI_COMM_WORLD);

            MPI_Send(&id, 1, MPI_INT, pid, MSG_TAG, MPI_COMM_WORLD);
        }

        /* Receive results from the workers */
        double local_determinant[m.row];
        for (int pid = 1, id = 0; pid < size; pid++)
        {
            MPI_Recv(&id, 1, MPI_INT, pid, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&local_determinant[id], 1, MPI_DOUBLE, pid, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        clock_gettime(CLOCK_MONOTONIC, &finish);
        long long diff = (finish.tv_sec - start.tv_sec) * 1000000000 + (finish.tv_nsec - start.tv_nsec);
        printf("Time elapsed was %lld ns \n", diff);

        fp = fopen("output.txt", "w");
        for (size_t i = 0; i < m.row; ++i)
        {
            double answer = local_determinant[i] / main_determinant;
            fprintf(fp, "%zu) %.2lf\n", i + 1, answer);
        }
        fclose(fp);

        free(m.data);
    }
    else
    { /* child process */

        struct Matrix m;

        MPI_Recv(&m.row, 1, MPI_INT, PARENT_ID, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&m.col, 1, MPI_INT, PARENT_ID, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        m.data = (double *)malloc(m.row * m.col * sizeof(double));
        MPI_Recv(&(m.data[0]), m.row * m.col, MPI_DOUBLE, PARENT_ID, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        m.activeCol = m.col - 1;

        int id = 0;
        MPI_Recv(&id, 1, MPI_INT, PARENT_ID, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        m.passedLine = (bool *)malloc(m.row * sizeof(bool));
        memset(m.passedLine, 0, m.row * sizeof(bool));
        double local_determinant = determinant(&m, 0, id); // CALCULATION

        MPI_Send(&id, 1, MPI_INT, PARENT_ID, MSG_TAG, MPI_COMM_WORLD);
        MPI_Send(&local_determinant, 1, MPI_DOUBLE, PARENT_ID, MSG_TAG, MPI_COMM_WORLD);

        /* Perhaps this is redundant */
        free(m.passedLine);
        free(m.data);
    }

    /* completes MPI */
    MPI_Finalize();
    return EXIT_SUCCESS;
}
