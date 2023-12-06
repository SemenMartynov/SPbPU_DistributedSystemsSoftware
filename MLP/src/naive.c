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

    for (int j = 0; j < n; j++) // столбцы исходной матрицы
    {
        int i_temp = 0;             // строка временной матрицы
        for (int i = 1; i < n; i++) // строки исходной матрицы
        {
            int j_temp = 0; // столбец временной матрицы
            // заполняем верхнюю половину временной матрицы
            for (int it = 0; it < j; it++)
            {
                temp[i_temp][j_temp++] = matrix[i][it];
            }
            // заполняем нижнюю половину временной матрицы
            // пропустив j-ый столбец
            for (int it = j + 1; it < n; it++)
            {
                temp[i_temp][j_temp++] = matrix[i][it];
            }
            i_temp++;
        }
        // временная матрица temp без 1-ой строки и j-го столбца
        det += sign * matrix[0][j] * determinant(temp, n - 1);
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