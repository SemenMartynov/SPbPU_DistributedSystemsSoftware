#include <stdio.h>
#include <mpi.h>
int main(int argc, char *argv[])
{
    int size1 = 0;
    int rank1 = 0;
    int size2 = 0;
    int rank2 = 0;
    int value = 0;

    MPI_Comm comm1;
    MPI_Status status;

    MPI_Init(&argc, &argv);                /* starts MPI */
    MPI_Comm_size(MPI_COMM_WORLD, &size1); /* get total process numer */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank1); /* get number of processes */

    printf("Всего в MPI_COMM_WORLD %d процесоов\n", size1);
    printf("Мой ранг в MPI_COMM_WORLD %d \n", rank1);

    MPI_Comm_split(MPI_COMM_WORLD, (rank1 % 2 == 0) ? 0 : 1, 0, &comm1);

    MPI_Comm_size(comm1, &size2);
    MPI_Comm_rank(comm1, &rank2);

    // 0-ой процесс нового коммуникатора, объединяющего все четные процессы
    if ((rank2 == 0) && (rank1 % 2))
        value = 777;
    
    // 0-ой процесс нового коммуникатора, объединяющего все нечетные процессы
    if ((rank2 == 0) && !(rank1 % 2))
        value = 666;

    MPI_Bcast(&value, 1, MPI_INT, 0, comm1);

    printf("Мой ранг в MPI_COMM_WORLD %d, Мой ранг в MPI_COMM_WORLD %d, broadcasted message %d", rank1, rank2, value);

    /* completes MPI */
    MPI_Finalize();
}