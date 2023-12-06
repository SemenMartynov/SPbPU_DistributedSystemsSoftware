# Многоязыковое параллельное программирование

## Постановка задачи

Задача: реализовать решение системы линейных уравнений (без ограничения на количество уравнений) методом Крамера.

В процессе работы необходимо получить и закрепить следующие компетенции:
- Разработка оптимизированного машинного кода;
- разработка векторизованного кода;
- разработка кода с выравниванием адресов;
- разработка собственного кода целевой платформы с использованием неявных средств многопоточности;
- разработка динамически подключаемых библиотек;
- кросскомпиляция кода с инструментальной платформы на целевую;
- многоязыковое программирование.

Для разработки подключаемой библиотеки использовать язык ​C11 (ISO/IEC 9899:2011). В качестве языка для внешней программы -- язык Rust (версия 1.72.0). Для неявной многопоточности использовать OpenMP не ниже версии 4.5.

## Тестовая программа

Тестовая программа, осуществляющая ввод/вывод данных и динамическую компоновку с библиотекой, содержащей функцию расчёта определителя.

```c
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

typedef double (*DeterminantFunc)(double **, int);

int main()
{
    int matrixSize;
    double **matrix;
    double *vector;

    printf("Введите количество уравнений: ");
    scanf("%d", &matrixSize);

    // Выделение памяти для матрицы системы
    matrix = (double **)malloc(matrixSize * sizeof(double *));
    for (int i = 0; i < matrixSize; i++)
    {
        matrix[i] = (double *)malloc(matrixSize * sizeof(double));
    }

    printf("Введите коэффициенты матрицы системы уравнений:\n");
    for (int i = 0; i < matrixSize; i++)
    {
        for (int j = 0; j < matrixSize; j++)
        {
            printf("Введите элемент [%d][%d]: ", i + 1, j + 1);
            scanf("%lf", &matrix[i][j]);
        }
    }

    // Выделение памяти для вектора свободных членов
    vector = (double *)malloc(matrixSize * sizeof(double));

    printf("Введите вектор свободных членов:\n");
    for (int i = 0; i < matrixSize; i++)
    {
        printf("Введите элемент [%d]: ", i + 1);
        scanf("%lf", &vector[i]);
    }

    // Загрузка динамической библиотеки
    void *library = dlopen("./libdeterminant.so", RTLD_LAZY);
    if (!library)
    {
        printf("Не удалось загрузить библиотеку\n");
        return 1;
    }

    // Получение указателя на функцию расчета определителя
    DeterminantFunc determinant = (DeterminantFunc)dlsym(library, "determinant");

    if (!determinant)
    {
        printf("Не удалось получить указатель на функцию\n");
        dlclose(library);
        return 1;
    }

    // Расчет определителя
    double det = determinant(matrix, matrixSize);

    // Проверка на нулевой определитель
    if (det == 0)
    {
        printf("Система уравнений не имеет единственного решения.\n");
        dlclose(library);
        return 1;
    }

    // Выделение памяти для решений
    double *root = (double *)malloc(matrixSize * sizeof(double));
    double **temp = (double **)malloc(matrixSize * sizeof(double *));
    for (int i = 0; i < matrixSize; i++)
    {
        temp[i] = (double *)malloc(matrixSize * sizeof(double));
    }

    for (int i = 0; i < matrixSize; i++)
    {
        // Создаем временную копию матрицы и заменяем i-ый столбец на вектор свободных членов
        for (int j = 0; j < matrixSize; j++)
        {
            for (int k = 0; k < matrixSize; k++)
            {
                temp[j][k] = matrix[j][k];
            }
        }

        for (int j = 0; j < matrixSize; j++)
        {
            temp[j][i] = vector[j];
        }

        // Расчет определителя для временной матрицы
        double det_temp = determinant(temp, matrixSize);

        // Решение по формуле xi = det(temp) / det(matrix)
        root[i] = det_temp / det;

        // Вывод решений
        printf("Решения системы:\n");
        printf("x%d = %.2f\n", i + 1, root[i]);
    }

    // Освобождение памяти и закрытие динамической библиотеки
    for (int i = 0; i < matrixSize; i++)
    {
        free(matrix[i]);
        free(temp[i]);
    }
    free(matrix);
    free(temp);
    free(vector);
    free(root);

    dlclose(library);

    return 0;
}
```

## Однопоточный невекторизованный код

### Исходный код

Код реализует алгоритм обработки данных без использования векторизации и неявных средств многопоточности.

```c
#include <stdio.h>
#include <stdlib.h>

double determinant(double **matrix, int n)
{
    if (n == 1)
        return matrix[0][0];

    if (n == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double det = 0;
    double **temp = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        temp[i] = (double *)malloc(n * sizeof(double));
    }
    int sign = 1;

    for (int k = 0; k < n; k++)
    {
        // Получаем временную матрицу, исключая первую строку и k-ый столбец
        int i_temp = 0;
        for (int i = 1; i < n; i++)
        {
            int j_temp = 0;
            for (int j = 0; j < n; j++)
            {
                if (j != k)
                {
                    temp[i_temp][j_temp++] = matrix[i][j];
                }
            }
            i_temp++;
        }

        det += sign * matrix[0][k] * determinant(temp, n - 1);
        sign = -sign;
    }

    // Освобождение памяти
    for (int i = 0; i < n; i++)
    {
        free(temp[i]);
    }
    free(temp);

    return det;
}
```

### Тестирование

**Тестирование не оптимизированного кода**

Для исследования, будем использовать общее время программы, и попробуем улучшить наиболее вычислительно-сложное место.
```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC
$ clang main.c -o main -std=c99 -ldl
$ time -p cat input10.txt | ./main
```

Для определения цикломатической сложности (McCabe's Cyclomatic Number) используем cccc (a code counter for C and C++).
```sh
$ cccc naive.c
```

Perf = 2.09
Comp = 11
P = Perf / Comp * 10^5 = 19000.000000000

**Тестирование кода с системными методами оптимизации**

```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops
$ clang main.c -o main -std=c99 -ldl
$ time -p cat input10.txt | ./main
```

Perf = 1.37  
Comp = 11  
P = Perf / Comp * 10^5 = 12454.545454545

## Однопоточный векторизованный код

### Исходный код

После долгих попыток подобрать решение с векторизацией (использовалась библиотека интрисиков `immintrin.h`), мы приходим к выводу, что невозможно в силу особенностей самого алгоритма.

В таком случае, попробуем поручить векторизацию компилятору.

### Тестирование

```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops -mavx2 -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize
naive.c:16:29: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        temp[i] = (double *)malloc(n * sizeof(double));
                            ^
naive.c:14:5: remark: loop not vectorized: instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
    for (int i = 0; i < n; i++)
    ^
naive.c:14:5: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
naive.c:27:13: remark: loop not vectorized: value that could not be identified as reduction is used outside the loop [-Rpass-analysis=loop-vectorize]
            for (int j = 0; j < n; j++)
            ^
naive.c:31:21: remark: loop not vectorized: cannot identify array bounds [-Rpass-analysis=loop-vectorize]
                    temp[i_temp][j_temp++] = matrix[i][j];
                    ^
naive.c:27:13: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
            for (int j = 0; j < n; j++)
            ^
naive.c:44:9: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        free(temp[i]);
        ^
naive.c:42:5: remark: loop not vectorized: instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
    for (int i = 0; i < n; i++)
    ^
naive.c:42:5: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
```

Тут мы явно указали использование AVX2 и попросили выводить всю отладочную информацию:
- -Rpass=loop-vectorize Отображать циклы, которые были успешно векторизованы
- -Rpass-missed=loop-vectorize Отображать циклы, которые не были успешно векторизованы
- -Rpass-analysis=loop-vectorize Отображать причины, которые не позволили выполнить векторизацию

С точки зрения производительности, это не имело заметного влияния.

## Однопоточный векторизованный код с использованием выравнивания адресов

Возложим задачу по выравниванию на компилятор

**Выравнивание по границе 32 Байта (2^5=32)**

Выравние адресов функций
```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops -mllvm -align-all-functions=5
```

Perf = 1.27  
Comp = 11  
P = Perf / Comp * 10^5 = 11545.454545455

Выравние всех блоков
```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops -mllvm -align-all-blocks=5
```

Perf = 1.28  
Comp = 11  
P = Perf / Comp * 10^5 = 11636.363636364

**Выравнивание по границе 64 Байта (2^6=64)**

Выравние адресов функций
```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops -mllvm -align-all-functions=6
```

Perf = 1.34  
Comp = 11  
P = Perf / Comp * 10^5 = 12181.818181818

Выравние всех блоков
```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops -mllvm -align-all-blocks=5
```

Perf = 1.37  
Comp = 11  
P = Perf / Comp * 10^5 = 12454.545454545

**Выравнивание по границе 128 Байта (2^7=128)**

Выравние адресов функций
```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops -mllvm -align-all-functions=7
```

Perf = 1.36  
Comp = 11  
P = Perf / Comp * 10^5 = 12363.636363636

Выравние всех блоков
```sh
$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops -mllvm -align-all-blocks=5
```

Perf = 1.71  
Comp = 11  
P = Perf / Comp * 10^5 = 15545.454545455

## Однопоточный векторизованный код с использованием выравнивания адресов

```sh
~/Development/SPbPU/PO_RVS/SPbPU_DistributedSystemsSoftware/MLP/src (main)
smart@thinkpad$ clang naive.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -funroll-loops -mavx2 -Rpass=loop-vectorize  -Rpass-analysis=loop-vectorize -mllvm -align-all-functions=5
naive.c:16:29: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        temp[i] = (double *)malloc(n * sizeof(double));
                            ^
naive.c:14:5: remark: loop not vectorized: instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
    for (int i = 0; i < n; i++)
    ^
naive.c:33:13: remark: vectorized loop (vectorization width: 4, interleaved count: 4) [-Rpass=loop-vectorize]
            for (int it = j + 1; it < n; it++)
            ^
naive.c:33:13: remark: vectorized loop (vectorization width: 4, interleaved count: 4) [-Rpass=loop-vectorize]
naive.c:27:13: remark: vectorized loop (vectorization width: 4, interleaved count: 2) [-Rpass=loop-vectorize]
            for (int it = 0; it < j; it++)
            ^
naive.c:27:13: remark: vectorized loop (vectorization width: 4, interleaved count: 4) [-Rpass=loop-vectorize]
naive.c:47:9: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        free(temp[i]);
        ^
naive.c:45:5: remark: loop not vectorized: instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
    for (int i = 0; i < n; i++)
    ^
```

Perf = 0.75  
Comp = 11  
P = Perf / Comp * 10^5 = 6818


## Многопоточный код с использованием векторизации, выравнивания адресов и неявных средств многопоточности

### Исходный код

Исходный код содержит неявниые средства многопоточности -- OpenMP.

```c
#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

double determinant(double **matrix, int n)
{
    if (n == 1)
        return matrix[0][0];

    if (n == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double det = 0;
    double **temp = (double **)malloc(n * sizeof(double *));
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        temp[i] = (double *)malloc(n * sizeof(double));
    }
    int sign = 1;

#pragma omp parallel for
    for (int k = 0; k < n; k++)
    {
        // Получаем временную матрицу, исключая первую строку и k-ый столбец
        int i_temp = 0;
#pragma omp parallel for
        for (int i = 1; i < n; i++)
        {
            int j_temp = 0;
#pragma omp parallel for
            for (int j = 0; j < n; j++)
            {
                if (j != k)
                {
                    temp[i_temp][j_temp++] = matrix[i][j];
                }
            }
            i_temp++;
        }

        det += sign * matrix[0][k] * determinant(temp, n - 1);
        sign = -sign;
    }

// Освобождение памяти
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        free(temp[i]);
    }
    free(temp);

    return det;
}
```

### Тестирование

```sh
$ clang openmp.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -fopenmp -funroll-loops -mllvm -align-all-blocks=5 -mavx2 -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize
openmp.c:19:29: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        temp[i] = (double *)malloc(n * sizeof(double));
                            ^
openmp.c:16:5: remark: loop not vectorized: instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
    #pragma omp parallel for
    ^
openmp.c:16:5: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
openmp.c:28:9: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        #pragma omp parallel for
        ^
openmp.c:23:5: remark: loop not vectorized: could not determine number of loop iterations [-Rpass-analysis=loop-vectorize]
    #pragma omp parallel for
    ^
openmp.c:23:5: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
openmp.c:32:13: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
            #pragma omp parallel for
            ^
openmp.c:28:9: remark: loop not vectorized: could not determine number of loop iterations [-Rpass-analysis=loop-vectorize]
        #pragma omp parallel for
        ^
openmp.c:28:9: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
openmp.c:37:21: remark: loop not vectorized: cannot identify array bounds [-Rpass-analysis=loop-vectorize]
                    temp[i_temp][j_temp++] = matrix[i][j];
                    ^
openmp.c:32:13: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
            #pragma omp parallel for
            ^
openmp.c:51:9: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        free(temp[i]);
        ^
openmp.c:48:5: remark: loop not vectorized: could not determine number of loop iterations [-Rpass-analysis=loop-vectorize]
    #pragma omp parallel for
    ^
openmp.c:48:5: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
```

Векторизация сталкивается с теми же проблемами -- не удаётся определить количество итераций. Это особенность алгоритма и рекурсивного подхода.

Perf = 4.74  
Comp = 11  
P = Perf / Comp * 10^5 = 43090.909090909

Использование неявных средств многопоточности дало худший результат из за слишком агрессивного распараллеливания. Кроме того, алгоритм начал выдавать не верные ответы.

Попробуем другую версию кода

```c
#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

double determinant(double **matrix, int n)
{
    if (n == 1)
        return matrix[0][0];

    if (n == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double det = 0;
    double **temp = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        temp[i] = (double *)malloc(n * sizeof(double));
    }
    int sign = 1;

#pragma omp parallel for
    for (int k = 0; k < n; k++)
    {
        // Получаем временную матрицу, исключая первую строку и k-ый столбец
        int i_temp = 0;
        for (int i = 1; i < n; i++)
        {
            int j_temp = 0;
            for (int j = 0; j < n; j++)
            {
                if (j != k)
                {
                    temp[i_temp][j_temp++] = matrix[i][j];
                }
            }
            i_temp++;
        }

        det += sign * matrix[0][k] * determinant(temp, n - 1);
        sign = -sign;
    }

// Освобождение памяти
    for (int i = 0; i < n; i++)
    {
        free(temp[i]);
    }
    free(temp);

    return det;
}
```

Компиляция со всеми ключами

```sh
$ clang openmp2.c -o libdeterminant.so -std=c99 -shared -fPIC -O3 -march=native -fopenmp -funroll-loops -mllvm -align-all-blocks=5 -mavx2 -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize
openmp2.c:18:29: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        temp[i] = (double *)malloc(n * sizeof(double));
                            ^
openmp2.c:16:5: remark: loop not vectorized: instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
    for (int i = 0; i < n; i++)
    ^
openmp2.c:16:5: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
openmp2.c:47:9: remark: loop not vectorized: call instruction cannot be vectorized [-Rpass-analysis=loop-vectorize]
        free(temp[i]);
        ^
openmp2.c:45:5: remark: loop not vectorized: could not determine number of loop iterations [-Rpass-analysis=loop-vectorize]
    for (int i = 0; i < n; i++)
    ^
openmp2.c:45:5: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
openmp2.c:30:13: remark: loop not vectorized: value that could not be identified as reduction is used outside the loop [-Rpass-analysis=loop-vectorize]
            for (int j = 0; j < n; j++)
            ^
openmp2.c:34:21: remark: loop not vectorized: cannot identify array bounds [-Rpass-analysis=loop-vectorize]
                    temp[i_temp][j_temp++] = matrix[i][j];
                    ^
openmp2.c:30:13: remark: loop not vectorized [-Rpass-missed=loop-vectorize]
            for (int j = 0; j < n; j++)
            ^
```

Perf = 0.56  
Comp = 11  
P = Perf / Comp * 10^5 = 5090.909090909

## Многоязычная программа

### Исходный код

Приложение на Rust, которое использует динамическую библиотеку

```Rs
extern crate libc;

use std::io;
use std::os::raw::{c_double, c_int};

#[link(name = "determinant", kind = "dylib")]
extern "C" {
    fn determinant(matrix: *mut *mut c_double, n: c_int) -> c_double;
}

fn main() {
    // Вводим размерность матрицы с клавиатуры
    println!("Введите размерность системы:");
    let mut dimension = String::new();
    io::stdin()
        .read_line(&mut dimension)
        .expect("Не удалось прочитать строку");
    let dimension: usize = dimension.trim().parse().expect("Некорректное число");

    // Создаем векторы для хранения коэффициентов матрицы и свободных членов
    let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; dimension]; dimension];
    let mut constants: Vec<f64> = vec![0.0; dimension];

    // Вводим коэффициенты матрицы и свободные члены с клавиатуры
    println!("Введите коэффициенты матрицы:");
    for i in 0..dimension {
        for j in 0..dimension {
            let mut coefficient = String::new();
            io::stdin()
                .read_line(&mut coefficient)
                .expect("Не удалось прочитать строку");
            let coefficient: f64 = coefficient.trim().parse().expect("Некорректное число");
            matrix[i][j] = coefficient;
        }
    }

    println!("Введите свободные члены:");
    for i in 0..dimension {
        let mut constant = String::new();
        io::stdin()
            .read_line(&mut constant)
            .expect("Не удалось прочитать строку");
        let constant: f64 = constant.trim().parse().expect("Некорректное число");
        constants[i] = constant;
    }

    // Создаем C-совместимую матрицу, чтобы передать ее C-функции
    let mut c_matrix: Vec<Vec<c_double>> = matrix
        .iter()
        .map(|row| row.iter().copied().map(|x| x as c_double).collect())
        .collect();
    let c_matrix_ptr: *mut *mut c_double = c_matrix
        .iter_mut()
        .map(|row| row.as_mut_ptr())
        .collect::<Vec<_>>()
        .as_mut_ptr();

        //println!("{:?}", matrix);

    // Вычисляем определитель с помощью C-функции
    let matrix_determinant = unsafe { determinant(c_matrix_ptr, dimension as c_int) };

    // Проверяем, что определитель не равен нулю (иначе система не имеет единственного решения)
    if matrix_determinant.abs() < 1e-10 {
        println!("Система уравнений не имеет единственного решения");
        return;
    }

    // Создаем вектор для хранения решений
    let mut solutions: Vec<f64> = vec![0.0; dimension];

    // Решаем систему уравнений методом Крамера
    for i in 0..dimension {
        let mut modified_matrix = matrix.clone();
        
        // Заменяем i-ый столбец на вектор свободных членов
        for j in 0..dimension {
            modified_matrix[j][i] = constants[j];
        }
        
        let c_matrix_i_ptr: *mut *mut c_double = modified_matrix
            .iter_mut()
            .map(|row| row.as_mut_ptr())
            .collect::<Vec<_>>()
            .as_mut_ptr();
        
        // Вычисляем определитель измененной матрицы
        let determinant_i = unsafe { determinant(c_matrix_i_ptr, dimension as c_int) };
        
        // Вычисляем решение и добавляем его в вектор решений
        solutions[i] = determinant_i / matrix_determinant;
    }

    // Выводим решение
    println!("Решение системы уравнений:");
    for (i, solution) in solutions.iter().enumerate() {
        println!("x{} = {}", i + 1, solution);
    }
}
```

### Тестирование

```sh
$ cat src/input10.txt | cargo run --package mlp --bin mlp
```

Perf = 1.64  
Comp = 11  
P = Perf / Comp * 10^5 = 14909.090909091

## Заключение

Использование системных методов оптимизации позволило повысить эффективность почти на 70%.

При выравнии, наилучший результат показало выравние на 32 байта. Дальнейшее увеличение только ухудшало этот результат.

Решение с использованием неявных средств многопоточности изначально дало худший результат из за слишком агрессивного распараллеливания. После пересмотра стратегии, мы получили самый лучший результат. Его можно улучшить путём изменения тестовой программы, но мы будем использовать это решение в многоязыковом решении.

Решение на Rust так же показало не высокую производительность. Вероятно это связано с тем, что для Rust требуются свои ключи компилятора для генерации эффективного кода.

Наиболее произодительным решением оказалось использование векторизации, выравнивания адресов и неявных средств многопоточности с демо-программой.

Обновление:
Наибольшая произовдительность у векторизованного и выровненного кода: 6'000 против 19'000 у самой наивной реализации.