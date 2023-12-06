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