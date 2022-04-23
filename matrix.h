#pragma once

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef int scalar_t;

typedef struct matrix {
    int num_cols_rows;
    scalar_t* matrix;
} matrix_t;

matrix_t* create_matrix(int size);
void delete_matrix(matrix_t* matrix);
void copy_matrix(matrix_t* src, matrix_t* dest);

scalar_t get_matrix_element(matrix_t* matrix, int row, int col);
void set_matrix_element(matrix_t* matrix, int row, int col, scalar_t element);

void add_matrix(matrix_t* a, matrix_t* b, matrix_t* c);
void subtract_matrix(matrix_t* a, matrix_t* b, matrix_t* c);

void print_matrix(matrix_t* matrix);