#include "matrix.h"

#define bump_free(x)

extern void* bump_pool;

static void* bump_malloc(size_t size)
{
    void* allocated_space = bump_pool;
    bump_pool = (void*)((uint8_t*)bump_pool + size);

    return allocated_space;
}

matrix_t* create_matrix(int rows, int cols)
{
    matrix_t* matrix = bump_malloc(sizeof(matrix_t));
    matrix->num_rows = rows;
    matrix->num_cols = cols;
    matrix->matrix = bump_malloc(sizeof(scalar_t) * rows * cols);

    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            set_matrix_element(matrix, y, x, 0);
        }
    }

    return matrix;
}

void reset_matrix(int rows, int cols, matrix_t* matrix)
{
    bump_free(matrix->matrix);

    matrix->num_rows = rows;
    matrix->num_cols = cols;
    matrix->matrix = bump_malloc(sizeof(scalar_t) * rows * cols);

    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            set_matrix_element(matrix, y, x, 0);
        }
    }
}

void delete_matrix(matrix_t* matrix)
{
    bump_free(matrix->matrix);
    bump_free(matrix);
}

scalar_t get_matrix_element(matrix_t* matrix, int row, int col)
{
    return matrix->matrix[(row * matrix->num_cols) + col];
}

void set_matrix_element(matrix_t* matrix, int row, int col, scalar_t element)
{
    matrix->matrix[(row * matrix->num_cols) + col] = element;
}

void copy_matrix(matrix_t* src, matrix_t* dest)
{
    for (int y = 0; y < src->num_rows; y++) {
        for (int x = 0; x < src->num_cols; x++) {
            set_matrix_element(dest, y, x, get_matrix_element(src, y, x));
        }
    }
}

void copy_to_smaller_matrix(matrix_t* src, matrix_t* dest)
{
    for (int y = 0; y < dest->num_rows; y++) {
        for (int x = 0; x < dest->num_cols; x++) {
            set_matrix_element(dest, y, x, get_matrix_element(src, y, x));
        }
    }
}

void add_matrix(matrix_t* a, matrix_t* b, matrix_t* c)
{
#ifdef DEBUG_TRANSFORMATIONS
    assert(!(a->num_rows != b->num_rows || a->num_cols != b->num_cols));
#endif
    for (int y = 0; y < a->num_rows; y++) {
        for (int x = 0; x < a->num_cols; x++) {
            scalar_t sum = get_matrix_element(a, y, x) + get_matrix_element(b, y, x);
            set_matrix_element(c, y, x, sum);
        }
    }
}

void subtract_matrix(matrix_t* a, matrix_t* b, matrix_t* c)
{
#ifdef DEBUG_TRANSFORMATIONS
    assert(!(a->num_rows != b->num_rows || a->num_cols != b->num_cols));
#endif
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