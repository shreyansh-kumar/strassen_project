#include "matrix.h"
#include "multiplication.h"
#include <stdlib.h>

static void fill_matrix(matrix_t* matrix, int rows, int cols, int seed);

int main(int argc, char* argv[])
{
    const int num_cols_rows = 1024;
    const int k = 16;
    matrix_t* a = create_matrix(num_cols_rows);
    matrix_t* b = create_matrix(num_cols_rows);
    matrix_t* res_a_b;

    fill_matrix(a, num_cols_rows, num_cols_rows, 1);
    fill_matrix(b, num_cols_rows, num_cols_rows, 2);

#if 0
    printf("a:\n");
    print_matrix(a);
    printf("\nb:\n");
    print_matrix(b);
#endif

    res_a_b = matrix_multiply(a, b, k);
#if 0
    printf("\na * b:\n");
    print_matrix(res_a_b);
#endif

    delete_matrix(a);
    delete_matrix(b);
    delete_matrix(res_a_b);

    return 0;
}

void fill_matrix(matrix_t* matrix, int rows, int cols, int seed)
{
    srand(seed);
    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            set_matrix_element(matrix, y, x, rand());
        }
    }
}