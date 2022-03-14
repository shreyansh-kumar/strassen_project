#include "multiplication.h"
#include "matrix.h"

int round_up_to_nearest_power_of_2(int num)
{
    num--;
    num |= num >> 1;
    num |= num >> 2;
    num |= num >> 4;
    num |= num >> 8;
    num |= num >> 16;
    num++;

    return num;
}

void copy_matrix_quadrants(
    matrix_t* src,
    matrix_t* src11,
    matrix_t* src12,
    matrix_t* src21,
    matrix_t* src22)
{
    int final_rows = src->num_rows;
    int final_cols = src->num_cols;
    if (final_rows >= 2)
        final_rows /= 2;
    if (final_cols >= 2)
        final_cols /= 2;
    for (int y = 0; y < final_rows; y++) {
        for (int x = 0; x < final_cols; x++) {
            set_matrix_element(src11, y, x, get_matrix_element(src, y, x));
            set_matrix_element(src12, y, x, get_matrix_element(src, y, x + final_cols));
            set_matrix_element(src21, y, x, get_matrix_element(src, y + final_rows, x));
            set_matrix_element(src22, y, x, get_matrix_element(src, y + final_rows, x + final_cols));
        }
    }
}

void copy_matrix_quadrants_out(
    matrix_t* src11,
    matrix_t* src12,
    matrix_t* src21,
    matrix_t* src22,
    matrix_t* dest)
{
    int final_rows = dest->num_rows;
    int final_cols = dest->num_cols;
    if (final_rows >= 2)
        final_rows /= 2;
    if (final_cols >= 2)
        final_cols /= 2;
    for (int y = 0; y < final_rows; y++) {
        for (int x = 0; x < final_cols; x++) {
            set_matrix_element(dest, y, x, get_matrix_element(src11, y, x));
            set_matrix_element(dest, y, x + final_cols, get_matrix_element(src12, y, x));
            set_matrix_element(dest, y + final_rows, x, get_matrix_element(src21, y, x));
            set_matrix_element(dest, y + final_rows, x + final_cols, get_matrix_element(src22, y, x));
        }
    }
}

matrix_t* matrix_multiply(matrix_t* a, matrix_t* b)
{
    int rounded_rows_a = round_up_to_nearest_power_of_2(a->num_rows);
    int rounded_cols_a = round_up_to_nearest_power_of_2(a->num_cols);
    int rounded_rows_b = round_up_to_nearest_power_of_2(b->num_rows);
    int rounded_cols_b = round_up_to_nearest_power_of_2(b->num_cols);

    matrix_t* padded_a = create_matrix(rounded_rows_a, rounded_cols_a);
    matrix_t* padded_b = create_matrix(rounded_rows_b, rounded_cols_b);

    copy_matrix(a, padded_a);
    copy_matrix(b, padded_b);

    matrix_t* padded_c = create_matrix(rounded_rows_a, rounded_cols_b);
    matrix_multiply_recursive_strassen(padded_a, padded_b, padded_c);

    matrix_t* c = create_matrix(a->num_rows, b->num_cols);
    copy_to_smaller_matrix(padded_c, c);

    delete_matrix(padded_a);
    delete_matrix(padded_b);
    delete_matrix(padded_c);

    return c;
}

void matrix_multiply_iterative(matrix_t* a, matrix_t* b, matrix_t* c)
{
    if (c->num_rows != a->num_rows || c->num_cols != b->num_cols)
        reset_matrix(a->num_rows, b->num_cols, c);

    for (int y = 0; y < c->num_rows; y++) {
        for (int x = 0; x < c->num_cols; x++) {
            scalar_t val = 0;
            for (int i = 0; i < a->num_cols; i++) {
                val += get_matrix_element(a, y, i) * get_matrix_element(b, i, x);
            }
            set_matrix_element(c, y, x, val);
        }
    }
}

#define SHOULD_DIVIDE_OR_SHOULD_ITERATE(x) \
    if (x > 1)                             \
        x /= 2;                            \
    else                                   \
        multiply = matrix_multiply_iterative;

void matrix_multiply_recursive_strassen(matrix_t* a, matrix_t* b, matrix_t* c)
{
    if (c->num_rows != a->num_rows || c->num_cols != b->num_cols)
        reset_matrix(a->num_rows, b->num_cols, c);

    if (
        (a->num_cols <= 2 && a->num_rows <= 2) || (b->num_cols <= 2 && b->num_rows <= 2)) {
        matrix_multiply_iterative(a, b, c);
        return;
    }
    register void (*multiply)(matrix_t*, matrix_t*, matrix_t*) = matrix_multiply_recursive_strassen;

    int half_rows_a = a->num_rows;
    int half_cols_a = a->num_cols;
    int half_rows_b = b->num_rows;
    int half_cols_b = b->num_cols;
    int half_rows_c = c->num_rows;
    int half_cols_c = c->num_cols;
    SHOULD_DIVIDE_OR_SHOULD_ITERATE(half_rows_a)
    SHOULD_DIVIDE_OR_SHOULD_ITERATE(half_cols_a)
    SHOULD_DIVIDE_OR_SHOULD_ITERATE(half_rows_b)
    SHOULD_DIVIDE_OR_SHOULD_ITERATE(half_cols_b)
    SHOULD_DIVIDE_OR_SHOULD_ITERATE(half_rows_c)
    SHOULD_DIVIDE_OR_SHOULD_ITERATE(half_cols_c)

    matrix_t* a11 = create_matrix(half_rows_a, half_cols_a);
    matrix_t* a12 = create_matrix(half_rows_a, half_cols_a);
    matrix_t* a21 = create_matrix(half_rows_a, half_cols_a);
    matrix_t* a22 = create_matrix(half_rows_a, half_cols_a);

    matrix_t* b11 = create_matrix(half_rows_b, half_cols_b);
    matrix_t* b12 = create_matrix(half_rows_b, half_cols_b);
    matrix_t* b21 = create_matrix(half_rows_b, half_cols_b);
    matrix_t* b22 = create_matrix(half_rows_b, half_cols_b);

    matrix_t* c11 = create_matrix(half_rows_c, half_cols_c);
    matrix_t* c12 = create_matrix(half_rows_c, half_cols_c);
    matrix_t* c21 = create_matrix(half_rows_c, half_cols_c);
    matrix_t* c22 = create_matrix(half_rows_c, half_cols_c);

    copy_matrix_quadrants(a, a11, a12, a21, a22);
    copy_matrix_quadrants(b, b11, b12, b21, b22);

#ifdef DEBUG_STRASSEN_RECURSION
    multiply = matrix_multiply_iterative;
#endif

    /*
     *
     * Strap yourself in, this is gonna get really hairy.
     *
     * Implementation lifted from:
     * S. Huss-Lederman, E. M. Jacobson, J. R. Johnson, A. Tsao, and T. Turnbull. Strassen's
     * algorithm for matrix multiplication: Modeling, analysis, and implementation. Technical
     * report, Center for Computing Sciences, 1996. Technical Report CCS-TR-96-147.
     *
     */

    // 1
    matrix_t* r1 = create_matrix(a11->num_rows, fmax(b11->num_rows, b11->num_cols));
#ifdef DEBUG_TRANSFORMATIONS
    printf("1 PRIOR <----------------------->\n");
    print_matrix(a11);
    printf("a11 has size %d-%d\n", a11->num_rows, a11->num_cols);
    print_matrix(a21);
    printf("a21 has size %d-%d\n", a21->num_rows, a21->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("1 PRIOR <----------------------->\n");
#endif
    subtract_matrix(a11, a21, r1);
#ifdef DEBUG_TRANSFORMATIONS
    printf("1 <----------------------->\n");
    print_matrix(a11);
    printf("a11 has size %d-%d\n", a11->num_rows, a11->num_cols);
    print_matrix(a21);
    printf("a21 has size %d-%d\n", a21->num_rows, a21->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("1 <----------------------->\n");
#endif

    // 2
    matrix_t* r2 = create_matrix(b22->num_rows, b22->num_cols);
#ifdef DEBUG_TRANSFORMATIONS
    printf("2 PRIOR <----------------------->\n");
    print_matrix(b22);
    printf("b22 has size %d-%d\n", b22->num_rows, b22->num_cols);
    print_matrix(b12);
    printf("b12 has size %d-%d\n", b12->num_rows, b12->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    printf("2 PRIOR <----------------------->\n");
#endif
    subtract_matrix(b22, b12, r2);
#ifdef DEBUG_TRANSFORMATIONS
    printf("2 <----------------------->\n");
    print_matrix(b22);
    printf("b22 has size %d-%d\n", b22->num_rows, b22->num_cols);
    print_matrix(b12);
    printf("b12 has size %d-%d\n", b12->num_rows, b12->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    printf("2 <----------------------->\n");
#endif

    // 3
    matrix_t* r3 = create_matrix(r1->num_rows, r2->num_cols);
#ifdef DEBUG_TRANSFORMATIONS
    printf("3 PRIOR <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r3);
    printf("r3 has size %d-%d\n", r3->num_rows, r3->num_cols);
    printf("3 PRIOR <----------------------->\n");
#endif
    multiply(r1, r2, r3);
#ifdef DEBUG_TRANSFORMATIONS
    printf("3 <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r3);
    printf("r3 has size %d-%d\n", r3->num_rows, r3->num_cols);
    printf("3 <----------------------->\n");
#endif

    // 4
#ifdef DEBUG_TRANSFORMATIONS
    printf("4 PRIOR <----------------------->\n");
    print_matrix(a21);
    printf("a21 has size %d-%d\n", a21->num_rows, a21->num_cols);
    print_matrix(a22);
    printf("a22 has size %d-%d\n", a22->num_rows, a22->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("4 PRIOR <----------------------->\n");
#endif
    add_matrix(a21, a22, r1);
#ifdef DEBUG_TRANSFORMATIONS
    printf("4 <----------------------->\n");
    print_matrix(a21);
    printf("a21 has size %d-%d\n", a21->num_rows, a21->num_cols);
    print_matrix(a22);
    printf("a22 has size %d-%d\n", a22->num_rows, a22->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("4 <----------------------->\n");
#endif

    // 5
#ifdef DEBUG_TRANSFORMATIONS
    printf("5 PRIOR <----------------------->\n");
    print_matrix(b12);
    printf("b12 has size %d-%d\n", b12->num_rows, b12->num_cols);
    print_matrix(b11);
    printf("b11 has size %d-%d\n", b11->num_rows, b11->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    printf("5 PRIOR <----------------------->\n");
#endif
    subtract_matrix(b12, b11, r2);
#ifdef DEBUG_TRANSFORMATIONS
    printf("5 <----------------------->\n");
    print_matrix(b12);
    printf("b12 has size %d-%d\n", b12->num_rows, b12->num_cols);
    print_matrix(b11);
    printf("b11 has size %d-%d\n", b11->num_rows, b11->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    printf("5 <----------------------->\n");
#endif

    // 6
    matrix_t* r4 = create_matrix(r1->num_rows, r2->num_cols);
#ifdef DEBUG_TRANSFORMATIONS
    printf("6 PRIOR <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r4);
    printf("r4 has size %d-%d\n", r4->num_rows, r4->num_cols);
    printf("6 PRIOR <----------------------->\n");
#endif
    multiply(r1, r2, r4);
#ifdef DEBUG_TRANSFORMATIONS
    printf("6 <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r4);
    printf("r4 has size %d-%d\n", r4->num_rows, r4->num_cols);
    printf("6 <----------------------->\n");
#endif

    // 7
#ifdef DEBUG_TRANSFORMATIONS
    printf("7 PRIOR <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(a11);
    printf("a11 has size %d-%d\n", a11->num_rows, a11->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("7 PRIOR <----------------------->\n");
#endif
    subtract_matrix(r1, a11, r1);
#ifdef DEBUG_TRANSFORMATIONS
    printf("7 <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(a11);
    printf("a11 has size %d-%d\n", a11->num_rows, a11->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("7 <----------------------->\n");
#endif

    // 8
#ifdef DEBUG_TRANSFORMATIONS
    printf("8 PRIOR <----------------------->\n");
    print_matrix(b22);
    printf("b22 has size %d-%d\n", b22->num_rows, b22->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    printf("8 PRIOR <----------------------->\n");
#endif
    subtract_matrix(b22, r2, r2);
#ifdef DEBUG_TRANSFORMATIONS
    printf("8 <----------------------->\n");
    print_matrix(b22);
    printf("b22 has size %d-%d\n", b22->num_rows, b22->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    printf("8 <----------------------->\n");
#endif

    // 9
    matrix_t* r5 = create_matrix(r1->num_rows, r2->num_cols);
#ifdef DEBUG_TRANSFORMATIONS
    printf("9 PRIOR <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    printf("9 PRIOR <----------------------->\n");
#endif
    multiply(r1, r2, r5);
#ifdef DEBUG_TRANSFORMATIONS
    printf("9 <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    printf("9 <----------------------->\n");
#endif

    // 10
#ifdef DEBUG_TRANSFORMATIONS
    printf("10 PRIOR <----------------------->\n");
    print_matrix(a12);
    printf("a12 has size %d-%d\n", a12->num_rows, a12->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("10 PRIOR <----------------------->\n");
#endif
    subtract_matrix(a12, r1, r1);
#ifdef DEBUG_TRANSFORMATIONS
    printf("10 <----------------------->\n");
    print_matrix(a12);
    printf("a12 has size %d-%d\n", a12->num_rows, a12->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("10 <----------------------->\n");
#endif

    // 11
#ifdef DEBUG_TRANSFORMATIONS
    printf("11 PRIOR <----------------------->\n");
    print_matrix(b21);
    printf("b21 has size %d-%d\n", b21->num_rows, b21->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    printf("11 PRIOR <----------------------->\n");
#endif
    subtract_matrix(b21, r2, r2);
#ifdef DEBUG_TRANSFORMATIONS
    printf("11 <----------------------->\n");
    print_matrix(b21);
    printf("b21 has size %d-%d\n", b21->num_rows, b21->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    printf("11 <----------------------->\n");
#endif

    // 12
    matrix_t* r6 = create_matrix(a22->num_rows, b22->num_cols);
#ifdef DEBUG_TRANSFORMATIONS
    printf("12 PRIOR <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(b22);
    printf("b22 has size %d-%d\n", b22->num_rows, b22->num_cols);
    print_matrix(r6);
    printf("r6 has size %d-%d\n", r6->num_rows, r6->num_cols);
    printf("12 PRIOR <----------------------->\n");
#endif
    multiply(r1, b22, r6);
#ifdef DEBUG_TRANSFORMATIONS
    printf("12 <----------------------->\n");
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(b22);
    printf("b22 has size %d-%d\n", b22->num_rows, b22->num_cols);
    print_matrix(r6);
    printf("r6 has size %d-%d\n", r6->num_rows, r6->num_cols);
    printf("12 <----------------------->\n");
#endif

    // 13
#ifdef DEBUG_TRANSFORMATIONS
    printf("13 PRIOR <----------------------->\n");
    print_matrix(r6);
    printf("r6 has size %d-%d\n", r6->num_rows, r6->num_cols);
    print_matrix(r4);
    printf("r4 has size %d-%d\n", r4->num_rows, r4->num_cols);
    print_matrix(r6);
    printf("r6 has size %d-%d\n", r6->num_rows, r6->num_cols);
    printf("13 PRIOR <----------------------->\n");
#endif
    add_matrix(r6, r4, r6);
#ifdef DEBUG_TRANSFORMATIONS
    printf("13 <----------------------->\n");
    print_matrix(r6);
    printf("r6 has size %d-%d\n", r6->num_rows, r6->num_cols);
    print_matrix(r4);
    printf("r4 has size %d-%d\n", r4->num_rows, r4->num_cols);
    print_matrix(r6);
    printf("r6 has size %d-%d\n", r6->num_rows, r6->num_cols);
    printf("13 <----------------------->\n");
#endif

    // 14
#ifdef DEBUG_TRANSFORMATIONS
    printf("14 PRIOR <----------------------->\n");
    print_matrix(a11);
    printf("a11 has size %d-%d\n", a11->num_rows, a11->num_cols);
    print_matrix(b11);
    printf("b11 has size %d-%d\n", b11->num_rows, b11->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("14 PRIOR <----------------------->\n");
#endif
    multiply(a11, b11, r1);
#ifdef DEBUG_TRANSFORMATIONS
    printf("14 <----------------------->\n");
    print_matrix(a11);
    printf("a11 has size %d-%d\n", a11->num_rows, a11->num_cols);
    print_matrix(b11);
    printf("b11 has size %d-%d\n", b11->num_rows, b11->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("14 <----------------------->\n");
#endif

    // 15
#ifdef DEBUG_TRANSFORMATIONS
    printf("15 PRIOR <----------------------->\n");
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    printf("15 PRIOR <----------------------->\n");
#endif
    add_matrix(r5, r1, r5);
#ifdef DEBUG_TRANSFORMATIONS
    printf("15 <----------------------->\n");
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    printf("15 <----------------------->\n");
#endif

    // 16
#ifdef DEBUG_TRANSFORMATIONS
    printf("16 PRIOR <----------------------->\n");
    print_matrix(r6);
    printf("r6 has size %d-%d\n", r6->num_rows, r6->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(c12);
    printf("c12 has size %d-%d\n", c12->num_rows, c12->num_cols);
    printf("16 PRIOR <----------------------->\n");
#endif
    add_matrix(r6, r5, r6);
    add_matrix(r6, c12, r6);
#ifdef DEBUG_TRANSFORMATIONS
    printf("16 <----------------------->\n");
    print_matrix(r6);
    printf("r6 has size %d-%d\n", r6->num_rows, r6->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(c12);
    printf("c12 has size %d-%d\n", c12->num_rows, c12->num_cols);
    printf("16 <----------------------->\n");
#endif

    // 17
#ifdef DEBUG_TRANSFORMATIONS
    printf("17 PRIOR <----------------------->\n");
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(r3);
    printf("r3 has size %d-%d\n", r3->num_rows, r3->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    printf("17 PRIOR <----------------------->\n");
#endif
    add_matrix(r5, r3, r5);
#ifdef DEBUG_TRANSFORMATIONS
    printf("17 <----------------------->\n");
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(r3);
    printf("r3 has size %d-%d\n", r3->num_rows, r3->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    printf("17 <----------------------->\n");
#endif

    // 18
#ifdef DEBUG_TRANSFORMATIONS
    printf("18 PRIOR <----------------------->\n");
    print_matrix(r4);
    printf("r4 has size %d-%d\n", r4->num_rows, r4->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(c22);
    printf("c22 has size %d-%d\n", c22->num_rows, c22->num_cols);
    printf("18 PRIOR <----------------------->\n");
#endif
    add_matrix(r4, r5, r4);
    add_matrix(r4, c22, r4);
#ifdef DEBUG_TRANSFORMATIONS
    printf("18 <----------------------->\n");
    print_matrix(r4);
    printf("r4 has size %d-%d\n", r4->num_rows, r4->num_cols);
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(c22);
    printf("c22 has size %d-%d\n", c22->num_rows, c22->num_cols);
    printf("18 <----------------------->\n");
#endif

    // 19
#ifdef DEBUG_TRANSFORMATIONS
    printf("19 PRIOR <----------------------->\n");
    print_matrix(a12);
    printf("a12 has size %d-%d\n", a12->num_rows, a12->num_cols);
    print_matrix(b21);
    printf("b21 has size %d-%d\n", b21->num_rows, b21->num_cols);
    print_matrix(r3);
    printf("r3 has size %d-%d\n", r3->num_rows, r3->num_cols);
    printf("19 PRIOR <----------------------->\n");
#endif
    multiply(a12, b21, r3);
#ifdef DEBUG_TRANSFORMATIONS
    printf("19 <----------------------->\n");
    print_matrix(a12);
    printf("a12 has size %d-%d\n", a12->num_rows, a12->num_cols);
    print_matrix(b21);
    printf("b21 has size %d-%d\n", b21->num_rows, b21->num_cols);
    print_matrix(r3);
    printf("r3 has size %d-%d\n", r3->num_rows, r3->num_cols);
    printf("19 <----------------------->\n");
#endif

    // 20
#ifdef DEBUG_TRANSFORMATIONS
    printf("20 PRIOR <----------------------->\n");
    print_matrix(r3);
    printf("r3 has size %d-%d\n", r3->num_rows, r3->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(c11);
    printf("c11 has size %d-%d\n", c11->num_rows, c11->num_cols);
    printf("20 PRIOR <----------------------->\n");
#endif
    add_matrix(r3, r1, r3);
    add_matrix(r3, c11, r3);
#ifdef DEBUG_TRANSFORMATIONS
    printf("20 <----------------------->\n");
    print_matrix(r3);
    printf("r3 has size %d-%d\n", r3->num_rows, r3->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(c11);
    printf("c11 has size %d-%d\n", c11->num_rows, c11->num_cols);
    printf("20 <----------------------->\n");
#endif

    // 21
#ifdef DEBUG_TRANSFORMATIONS
    printf("21 PRIOR <----------------------->\n");
    print_matrix(a22);
    printf("a22 has size %d-%d\n", a22->num_rows, a22->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("21 PRIOR <----------------------->\n");
#endif
    multiply(a22, r2, r1);
#ifdef DEBUG_TRANSFORMATIONS
    printf("21 <----------------------->\n");
    print_matrix(a22);
    printf("a22 has size %d-%d\n", a22->num_rows, a22->num_cols);
    print_matrix(r2);
    printf("r2 has size %d-%d\n", r2->num_rows, r2->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    printf("21 <----------------------->\n");
#endif

    // 22
#ifdef DEBUG_TRANSFORMATIONS
    printf("22 PRIOR <----------------------->\n");
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(c21);
    printf("c21 has size %d-%d\n", c21->num_rows, c21->num_cols);
    printf("22 PRIOR <----------------------->\n");
#endif
    add_matrix(r5, r1, r5);
    add_matrix(r5, c21, r5);
#ifdef DEBUG_TRANSFORMATIONS
    printf("22 <----------------------->\n");
    print_matrix(r5);
    printf("r5 has size %d-%d\n", r5->num_rows, r5->num_cols);
    print_matrix(r1);
    printf("r1 has size %d-%d\n", r1->num_rows, r1->num_cols);
    print_matrix(c21);
    printf("c21 has size %d-%d\n", c21->num_rows, c21->num_cols);
    printf("22 <----------------------->\n");
#endif

    copy_matrix_quadrants_out(r3, r6, r5, r4, c);

    delete_matrix(r1);
    delete_matrix(r2);
    delete_matrix(r3);
    delete_matrix(r4);
    delete_matrix(r5);
    delete_matrix(r6);

    delete_matrix(a11);
    delete_matrix(a12);
    delete_matrix(a21);
    delete_matrix(a22);

    delete_matrix(b11);
    delete_matrix(b12);
    delete_matrix(b21);
    delete_matrix(b22);

    delete_matrix(c11);
    delete_matrix(c12);
    delete_matrix(c21);
    delete_matrix(c22);
}