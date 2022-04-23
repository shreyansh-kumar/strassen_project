#include "matrix.h"
#include "multiplication.h"

static void fill_matrix(matrix_t* matrix, int rows, int cols);

int main(int argc, char* argv[])
{
    const int k = 16;
    matrix_t* a = create_matrix(1024);
    matrix_t* b = create_matrix(1024);
    matrix_t* res_a_b;

    fill_matrix(a, 8, 8);
    fill_matrix(b, 8, 4);

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

void fill_matrix(matrix_t* matrix, int rows, int cols)
{
    scalar_t val = 0;
    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            set_matrix_element(matrix, y, x, val);
            val++;
            if (val > 7) {
                val = 0;
            }
        }
    }
}