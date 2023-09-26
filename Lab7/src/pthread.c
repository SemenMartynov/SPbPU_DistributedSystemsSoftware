#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <float.h>

#include <string.h>
#include <pthread.h>

struct Matrix
{
    double *data;     // Contents of the matrix
    size_t row;       // Number of lines (equations)
    size_t col;       // Number of columns (including answers)
    size_t activeCol; // Number of columns
    bool *passedLine; // Excluded lines
};

struct Params
{
    size_t id;
    double local_determinant;
    struct Matrix *m;
};

double determinant(struct Matrix *m, size_t colIt, size_t colSub /*col to be substituted*/)
{
    if (colIt == m->activeCol)
        return 1;

    double sum = 0;
    int sign = 1;
    for (size_t i = 0; i < m->row; ++i)
    {
        // Recursive way to reduce matrix dimension
        if (!m->passedLine[i])
        {
            m->passedLine[i] = true;
            // When we get to the column that needs to be replaced,
            // we take the values ​​from the very last column with the answers.
            size_t shift = colIt == colSub ? m->activeCol : colIt;
            sum += sign * m->data[i * m->col + shift] * determinant(m, colIt + 1, colSub);
            m->passedLine[i] = false;
            sign *= -1;
        }
    }
    return sum;
}

void *_determinant(void *params)
{
    struct Params *p = (struct Params *)params;
    p->m->passedLine = (bool *)malloc(p->m->row * sizeof(bool));
    memset(p->m->passedLine, 0, p->m->row * sizeof(bool));
    p->local_determinant = determinant(p->m, 0, p->id);
    free(p->m->passedLine);
    pthread_exit(NULL);
}

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

    pthread_t threads[m.row];
    struct Params params[m.row];
    for (size_t i = 0; i != m.row; ++i)
    {
        params[i].id = i; // index
        params[i].m = &m; // matrix
    }
    // Run threads
    for (size_t i = 0; i != m.row; ++i)
    {
        pthread_create(&threads[i], NULL, _determinant, &params[i]);
    }
    // Waiting for threads to complete
    for (size_t i = 0; i != m.row; ++i)
    {
        pthread_join(threads[i], NULL);
    }

    // Getting roots and recording results
    fp = fopen("output.txt", "w");
    for (size_t i = 0; i != m.row; ++i)
    {
        double answer = params[i].local_determinant / main_determinant;
        fprintf(fp, "%zu) %.2lf\n", params[i].id + 1, answer);
    }
    free(m.data);
    fclose(fp);
    return EXIT_SUCCESS;
}