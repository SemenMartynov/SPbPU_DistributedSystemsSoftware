#include <stdbool.h>
#include "matrix.h"

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