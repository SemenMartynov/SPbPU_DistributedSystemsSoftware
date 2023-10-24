#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <float.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include <omp.h>

#include "matrix.h"

int main()
{
    struct Matrix m;
    // Read input from file
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
        return EXIT_SUCCESS;
    }

    double local_determinant[m.row];
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
#pragma omp parallel for firstprivate(m)
    for (size_t i = 0; i < m.row; ++i)
    {
        m.passedLine = (bool *)malloc(m.row * sizeof(bool)); // EXPENSIVE!
        memset(m.passedLine, 0, m.row * sizeof(bool));
        local_determinant[i] = determinant(&m, 0, i);
        free(m.passedLine);
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
    return EXIT_SUCCESS;
}