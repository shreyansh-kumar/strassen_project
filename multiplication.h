#pragma once

#include "matrix.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum MULT_FUNC {
    SAM,
    BAM
};

matrix_t* matrix_multiply(matrix_t* a, matrix_t* b, enum MULT_FUNC function, int k);
void matrix_multiply_iterative(matrix_t* a, matrix_t* b, matrix_t* c);
void matrix_multiply_recursive_strassen(matrix_t* a, matrix_t* b, matrix_t* c, int k);