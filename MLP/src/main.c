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
    double *x = (double *)malloc(matrixSize * sizeof(double));
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
        x[i] = det_temp / det;
    }

    // Вывод решений
    printf("Решения системы:\n");
    for (int i = 0; i < matrixSize; i++)
    {
        printf("x%d = %.2f\n", i + 1, x[i]);
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
    free(x);

    dlclose(library);

    return 0;
}