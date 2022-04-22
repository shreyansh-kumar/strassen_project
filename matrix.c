#include "matrix.h"

matrix_t* create_matrix(int rows, int cols)
{
    matrix_t* matrix = malloc(sizeof(matrix_t));
    matrix->num_rows = rows;
    matrix->num_cols = cols;
    matrix->matrix = malloc(sizeof(scalar_t) * rows * cols);

    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            set_matrix_element(matrix, y, x, 0);
        }
    }

    return matrix;
}

void delete_matrix(matrix_t* matrix)
{
    free(matrix->matrix);
    free(matrix);
}

void copy_matrix(matrix_t* src, matrix_t* dest)
{
    for (int y = 0; y < src->num_rows; y++) {
        for (int x = 0; x < src->num_cols; x++) {
            set_matrix_element(dest, y, x, get_matrix_element(src, y, x));
        }
    }
}

scalar_t get_matrix_element(matrix_t* matrix, int row, int col)
{
    return matrix->matrix[(row * matrix->num_cols) + col];
}

void set_matrix_element(matrix_t* matrix, int row, int col, scalar_t element)
{
    matrix->matrix[(row * matrix->num_cols) + col] = element;
}

void add_matrix(matrix_t* a, matrix_t* b, matrix_t* c)
{
    for (int y = 0; y < a->num_rows; y++) {
        for (int x = 0; x < a->num_cols; x++) {
            scalar_t sum = get_matrix_element(a, y, x) + get_matrix_element(b, y, x);
            set_matrix_element(c, y, x, sum);
        }
    }
}

void subtract_matrix(matrix_t* a, matrix_t* b, matrix_t* c)
{
    for (int y = 0; y < a->num_rows; y++) {
        for (int x = 0; x < a->num_cols; x++) {
            scalar_t difference = get_matrix_element(a, y, x) - get_matrix_element(b, y, x);
            set_matrix_element(c, y, x, difference);
        }
    }
}

void print_matrix(matrix_t* matrix)
{
    for (int y = 0; y < matrix->num_rows; y++) {
        for (int x = 0; x < matrix->num_cols; x++) {
            printf("%d ", get_matrix_element(matrix, y, x));
        }
        printf("\n");
    }
}