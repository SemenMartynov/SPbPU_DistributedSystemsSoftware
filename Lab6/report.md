# Лабораторная работа 6. Создание простого приложения в среде OpenMP

Цели и задачи:  
Научиться использовать OpenMP для разработки параллельных программ и распараллеливания уже написанных программ.

Вариант 8.  
Дана последовательность натуральных чисел {a[0]...a[n–1]}. Создать OpenMP-приложение для поиска суммы ∑a[i], где a[i] -- четные числа.

## Реализация алгоритма в последовательной программе

Реализация программы
```c
#include <stdio.h>
#include <stdbool.h>

int data[] = {14, 32, 36, 34, 95, 36, 57, 18, 19, 16};

static inline bool is_even(int x) { return (x & 1) == 0; }

int main()
{
    int sum = 0;

    size_t length = sizeof(data) / sizeof(data[0]);
    for (size_t i = 0; i != length; i++)
    {
        if (is_even(data[i]))
        {
            sum += data[i];
        }
    }

    printf("The sum of the even natural numbers is equal to %d.\n", sum);
    return 0;
}
```

Компиляция и запуск приложения
```sh
$ gcc -O2 sequential.c -o sequential
$ time ./sequential 
The sum of the even natural numbers is equal to 186.

real	0m0.003s
user	0m0.002s
sys	0m0.001s
```

## Распараллеливание программы при помощи директив OpenMP

Реализация программы
```c
#include <stdio.h>
#include <stdbool.h>

#include <omp.h>

int data[] = {14, 32, 36, 34, 95, 36, 57, 18, 19, 16};

static inline bool is_even(int x) { return (x & 1) == 0; }

int main()
{
    int sum = 0;

    size_t length = sizeof(data) / sizeof(data[0]);
#pragma omp parallel for reduction (+:sum)
    for (size_t i = 0; i != length; i++)
    {
        if (is_even(data[i]))
        {
            sum += data[i];
        }
    }

    printf("The sum of the even natural numbers is equal to %d.\n", sum);
    return 0;
}
```

Условие `reduction` позволяет производить безопасное глобальное вычисление. Приватная копия каждой перечисленной переменной инициализируется при входе в параллельную секцию в соответствии с указанным оператором (0 для оператора +). При выходе из параллельной секции из частично вычисленных значений вычисляется результирующее и передается в основной поток.

Задать количество потоков можно таким приёмом `omp_set_num_threads(omp_get_num_threads() < 4 ? omp_get_num_threads() : 4);`. Либо можно использовать переменную окружения `OMP_NUM_THREADS=4`. При этом функции внутри кода имеют бОльший приоритет.

Компиляция и запуск приложения
```sh
$ gcc -fopenmp -O2 parallel.c -o parallel
$ time OMP_NUM_THREADS=4 ./parallel
The sum of the even natural numbers is equal to 186.

real	0m0.004s
user	0m0.000s
sys	0m0.004s
```

## Выводы

В данной работе последовательная программа отработала быстрее параллельной, т.к. вычисления очень простые, а работа с потоками даёт накладные расходы.
Если проанализировать чем занималась программа, наверняка наибольшее время уйдёт не на вычисление, а на вывод результата.