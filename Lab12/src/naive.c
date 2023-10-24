#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include <float.h>
#include <stdint.h>
#include <time.h>

#include "matrix.h"

int main()
{
    struct Matrix m;
    // Read input from file
    FILE *fp = fopen("input.txt", "r");
    fscanf(fp, "%zu %zu", &m.row, &m.col);
    m.data = (double *)malloc(m.row * m.col * sizeof(double));
    m.passedLine = (bool *)malloc(m.row * sizeof(bool));
    for (size_t i = 0; i != m.row; ++i)
    {
        for (size_t j = 0; j != m.col; ++j)
        {
            fscanf(fp, "%lf", &m.data[i * m.col + j]);
        }
    }
    m.activeCol = m.col - 1;
    fclose(fp);

    double main_determinant = determinant(&m, 0, m.activeCol);
    if (fabs(main_determinant) < DBL_EPSILON)
    {
        printf("\nDetermining roots for a given system of equations is not possible.\n");
        return EXIT_SUCCESS;
    }

    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);

    double local_determinant[m.row];
    for (size_t sub = 0; sub < m.row; ++sub)
    {
        local_determinant[sub] = determinant(&m, 0, sub);
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);
    long long diff = (finish.tv_sec - start.tv_sec) * 1000000000 + (finish.tv_nsec - start.tv_nsec);
    printf("Time elapsed was %lld ns \n", diff);

    // Getting roots and recording results
    fp = fopen("output.txt", "w");
    for (size_t i = 0; i != m.row; ++i)
    {
        double answer = local_determinant[i] / main_determinant;
        fprintf(fp, "%zu) %.2lf\n", i + 1, answer);
    }
    fclose(fp);

    free(m.data);
    free(m.passedLine);
    return EXIT_SUCCESS;
}