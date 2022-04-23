#include "matrix.h"

extern size_t max_heap_space_used_bytes;
extern size_t current_heap_space_used_bytes;

matrix_t* create_matrix(int size)
{
    matrix_t* matrix = malloc(sizeof(matrix_t));
    matrix->num_cols_rows = size;
    matrix->matrix = malloc(sizeof(scalar_t) * size * size);

    current_heap_space_used_bytes += sizeof(matrix_t);
    current_heap_space_used_bytes += (sizeof(scalar_t) * size * size);

    if (current_heap_space_used_bytes > max_heap_space_used_bytes)
        max_heap_space_used_bytes = current_heap_space_used_bytes;

    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            set_matrix_element(matrix, y, x, 0);
        }
    }

    return matrix;
}

void delete_matrix(matrix_t* matrix)
{
    current_heap_space_used_bytes -= sizeof(matrix_t);
    current_heap_space_used_bytes -= (sizeof(scalar_t) * matrix->num_cols_rows * matrix->num_cols_rows);

    free(matrix->matrix);
    free(matrix);
}

void copy_matrix(matrix_t* src, matrix_t* dest)
{
    for (int y = 0; y < src->num_cols_rows; y++) {
        for (int x = 0; x < src->num_cols_rows; x++) {
            set_matrix_element(dest, y, x, get_matrix_element(src, y, x));
        }
    }
}

scalar_t get_matrix_element(matrix_t* matrix, int row, int col)
{
    return matrix->matrix[(row * matrix->num_cols_rows) + col];
}

void set_matrix_element(matrix_t* matrix, int row, int col, scalar_t element)
{
    matrix->matrix[(row * matrix->num_cols_rows) + col] = element;
}

void add_matrix(matrix_t* a, matrix_t* b, matrix_t* c)
{
    for (int y = 0; y < a->num_cols_rows; y++) {
        for (int x = 0; x < a->num_cols_rows; x++) {
            scalar_t sum = get_matrix_element(a, y, x) + get_matrix_element(b, y, x);
            set_matrix_element(c, y, x, sum);
        }
    }
}

void subtract_matrix(matrix_t* a, matrix_t* b, matrix_t* c)
{
    for (int y = 0; y < a->num_cols_rows; y++) {
        for (int x = 0; x < a->num_cols_rows; x++) {
            scalar_t difference = get_matrix_element(a, y, x) - get_matrix_element(b, y, x);
            set_matrix_element(c, y, x, difference);
        }
    }
}

void print_matrix(matrix_t* matrix)
{
    for (int y = 0; y < matrix->num_cols_rows; y++) {
        for (int x = 0; x < matrix->num_cols_rows; x++) {
            printf("%d ", get_matrix_element(matrix, y, x));
        }
        printf("\n");
    }
}