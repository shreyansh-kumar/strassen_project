#pragma once

#include "matrix.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void copy_matrix_quadrants(
    matrix_t* src,
    matrix_t* src11,
    matrix_t* src12,
    matrix_t* src21,
    matrix_t* src22);
void copy_matrix_quadrants_out(
    matrix_t* src11,
    matrix_t* src12,
    matrix_t* src21,
    matrix_t* src22,
    matrix_t* dest);
matrix_t* matrix_multiply(matrix_t* a, matrix_t* b);
void matrix_multiply_iterative(matrix_t* a, matrix_t* b, matrix_t* c);
void matrix_multiply_recursive_strassen(matrix_t* a, matrix_t* b, matrix_t* c);

int round_up_to_nearest_power_of_2(int num);