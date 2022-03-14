#include "matrix.h"
#include "multiplication.h"

static void* start_bump_pool;
void* bump_pool;

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

int main(int argc, char* argv[])
{
    start_bump_pool = malloc(32 * 1024 * 1024);
    bump_pool = start_bump_pool;

    matrix_t* a = create_matrix(8, 8);
    matrix_t* b = create_matrix(8, 4);
    matrix_t* res_a_b;

    fill_matrix(a, 8, 8);
    fill_matrix(b, 8, 4);

    printf("a:\n");
    print_matrix(a);
    printf("\nb:\n");
    print_matrix(b);

    res_a_b = matrix_multiply(a, b);
    printf("\na * b:\n");
    print_matrix(res_a_b);

    delete_matrix(a);
    delete_matrix(b);
    delete_matrix(res_a_b);

    free(start_bump_pool);

    return 0;
}