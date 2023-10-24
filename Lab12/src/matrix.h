#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>

struct Matrix
{
    double *data;     // Contents of the matrix
    size_t row;       // Number of lines (equations)
    size_t col;       // Number of columns (including answers)
    size_t activeCol; // Number of columns
    bool *passedLine; // Excluded lines
};

double determinant(struct Matrix *m, size_t colIt, size_t colSub);

#endif /* MATRIX_H */