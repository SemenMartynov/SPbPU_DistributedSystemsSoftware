#include <mpi.h>
#include <stdio.h>
#define leng 20 // length of WR-string

int main(int argc, char **argv)
{
    double t, t2;
    char WR[leng];

    int rank, size;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    sprintf(WR, "Hello from %d", rank); // формирование сообщения

    t = MPI_Wtime();                    // фиксация времени «начала посылки»,
                                        // локально для каждого процесса
    MPI_Send(WR, leng, MPI_CHAR, size - (rank + 1), rank, MPI_COMM_WORLD);          // (1)
    MPI_Recv(WR, leng, MPI_CHAR, rank, size - (rank + 1), MPI_COMM_WORLD, &status); // (2)
    t2 = MPI_Wtime();                   // фиксация времени «окончания приема»,
                                        // локально для каждого процесса
    printf("\n From processor %d\n WR=%s\n", rank, WR);                // вывод сообщения
    printf("\n From processor %d\n Time=%le\n", rank, (t2 - t) / 100); // вывод времени,
                                                    // затраченного на обмен данным процессором
    MPI_Finalize();
    return 0;
}