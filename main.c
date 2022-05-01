#include "matrix.h"
#include "multiplication.h"
#include <stdlib.h>

static void fill_matrix(matrix_t* matrix, int rows, int cols, int seed);

int main(int argc, char* argv[])
{
    const int num_cols_rows = 1024;
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

    // If you want to perform a SAM instead of a SAMk, set k to 1
    // If you want to perform a BAM, switch SAM to BAM and vice versa
    for (int k = 1; k <= 1024; k *= 2) {
        res_a_b = matrix_multiply(a, b, SAM, k);
        delete_matrix(res_a_b);
    }
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