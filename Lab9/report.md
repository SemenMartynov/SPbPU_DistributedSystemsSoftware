# Лабораторная работа 9. Создание простого приложения с помощью библиотеки MPI

Цели и задачи:  
Изучить различные способы двухточечного обмена в MPI,
производные типы данных и операции упаковки и распаковки данных.

Вариант 8. Аукцион.  

Известный предприниматель и филантроп Джон Смит непосильным трудом на ниве экономических преступлений нажил огромное состояние и умер. Имущество досталось его сыну Джону Смиту младшему, который, в отличие от отца, отличается отменным здоровьем и любовью к рулетке. Проиграв за полгода заводы, газеты и пароходы, приобретенные отцом, он для поправки дел решает распродать с аукциона вещи из отцовского особняка.

Аукцион проводится следующим образом: каждый из участников в тайне от остальных пишет свою цену на специальной карточке и отдает ее распорядителю. Просмотрев карточки, распорядитель объявляет победителя.

Написать программу, моделирующую проведение аукциона, используя метод
передачи информации «точка-точка».

## Приложение


Исходный код приложения
```c
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
```

Компиляция и запуск приложения на 6 ядрах

```sh
$ gcc -std=c99 -Wall -O2 openmp.c -o openmp -lmpi

$ mpiexec -n 6 openmp
The auction participant 4 bids 611
The auction participant 1 bids 799
The auction participant 3 bids 775
The auction participant 5 bids 190
The auction participant 2 bids 424
RReceived bid 799 from the auction participant 1
Received bid 424 from the auction participant 2
Received bid 775 from the auction participant 3
Received bid 611 from the auction participant 4
Received bid 190 from the auction participant 5
Max bid is 799, the participant 1 win!
```

## Заключение

В данной работе демонстрируется реализация алгоритма аукциона с использованием функций двухточечного обмена библиотеки MPI.  
В виду простоты программы, замеры времени не дадут интересных результатов, хотя достаточно очевидно, что вычислительно-сложные задачи выигрывают от такого распараллеливания.  
Что касается упаковки данных, тут она так же является избыточной, так как передаются случайно сгенерированные целочисленные значения.