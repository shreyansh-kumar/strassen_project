#include "multiplication.h"
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#define USEC_PER_SEC 1000000

size_t max_stack_location;
size_t min_stack_location;

static void copy_matrix_quadrants(
    matrix_t* src,
    matrix_t* src11,
    matrix_t* src12,
    matrix_t* src21,
    matrix_t* src22);

static void copy_matrix_quadrants_out(
    matrix_t* src11,
    matrix_t* src12,
    matrix_t* src21,
    matrix_t* src22,
    matrix_t* dest);

matrix_t* matrix_multiply(matrix_t* a, matrix_t* b, int k)
{
    // Sanity check that the matrices are compatible
    assert(a->num_cols_rows == b->num_cols_rows);

    matrix_t* c = create_matrix(a->num_cols_rows);

    // Find the start of the stack space used (on x86, the stack grows down)
    max_stack_location = (size_t)&a;
    min_stack_location = max_stack_location;

    /*
     * Ok so there's a caveat here, since we fork, this doesn't return
     * the matrix c to the parent process. Normally this would be a big
     * deal, but we aren't actually gonna use matrix c for anything, we're
     * purely just trying to measure the time and space. As such, we can just
     * safely discard c and just return a blank matrix back to the parent process
     */
    int child = fork();
    if (child == 0) {
        matrix_multiply_recursive_strassen(a, b, c, k);

        struct rusage usage;

        assert(getrusage(RUSAGE_SELF, &usage) == 0);

        size_t cpu_time_used_seconds_in_usec = (usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) * USEC_PER_SEC;
        size_t cpu_total_time_used_usec = cpu_time_used_seconds_in_usec + usage.ru_utime.tv_usec + usage.ru_stime.tv_usec;

        size_t total_stack_space_used_bytes = max_stack_location - min_stack_location;

        printf("c->num_cols_rows: %d -- k value: %d -- cpu_total_time_used_usec: %zu, total_stack_space_used_bytes: %zu\n",
            c->num_cols_rows,
            k,
            cpu_total_time_used_usec,
            total_stack_space_used_bytes);

        // Avoid any memory leaks
        delete_matrix(a);
        delete_matrix(b);
        delete_matrix(c);

        exit(EXIT_SUCCESS);
    } else {
        wait(NULL);
    }

    return c;
}

void matrix_multiply_iterative(matrix_t* a, matrix_t* b, matrix_t* c)
{
    // Figure out how much stack space we're currently using
    size_t curr_stack_location = (size_t)&c;
    if (curr_stack_location < min_stack_location)
        min_stack_location = curr_stack_location;

    for (int y = 0; y < c->num_cols_rows; y++) {
        for (int x = 0; x < c->num_cols_rows; x++) {
            scalar_t val = 0;
            for (int i = 0; i < a->num_cols_rows; i++) {
                val += get_matrix_element(a, y, i) * get_matrix_element(b, i, x);
            }
            set_matrix_element(c, y, x, val);
        }
    }
}

void matrix_multiply_recursive_strassen(matrix_t* a, matrix_t* b, matrix_t* c, int k)
{
    if (a->num_cols_rows <= k || b->num_cols_rows <= k) {
        matrix_multiply_iterative(a, b, c);
        return;
    }

    int matrix_quad_size = a->num_cols_rows / 2;

    matrix_t* a11 = create_matrix(matrix_quad_size);
    matrix_t* a12 = create_matrix(matrix_quad_size);
    matrix_t* a21 = create_matrix(matrix_quad_size);
    matrix_t* a22 = create_matrix(matrix_quad_size);

    matrix_t* b11 = create_matrix(matrix_quad_size);
    matrix_t* b12 = create_matrix(matrix_quad_size);
    matrix_t* b21 = create_matrix(matrix_quad_size);
    matrix_t* b22 = create_matrix(matrix_quad_size);

    matrix_t* c11 = create_matrix(matrix_quad_size);
    matrix_t* c12 = create_matrix(matrix_quad_size);
    matrix_t* c21 = create_matrix(matrix_quad_size);
    matrix_t* c22 = create_matrix(matrix_quad_size);

    copy_matrix_quadrants(a, a11, a12, a21, a22);
    copy_matrix_quadrants(b, b11, b12, b21, b22);

    matrix_t* s1 = create_matrix(matrix_quad_size);
    subtract_matrix(b12, b22, s1);

    matrix_t* s2 = create_matrix(matrix_quad_size);
    add_matrix(a11, a12, s2);

    matrix_t* s3 = create_matrix(matrix_quad_size);
    add_matrix(a21, a22, s3);

    matrix_t* s4 = create_matrix(matrix_quad_size);
    subtract_matrix(b21, b11, s4);

    matrix_t* s5 = create_matrix(matrix_quad_size);
    add_matrix(a11, a22, s5);

    matrix_t* s6 = create_matrix(matrix_quad_size);
    add_matrix(b11, b22, s6);

    matrix_t* s7 = create_matrix(matrix_quad_size);
    subtract_matrix(a12, a22, s7);

    matrix_t* s8 = create_matrix(matrix_quad_size);
    add_matrix(b21, b22, s8);

    matrix_t* s9 = create_matrix(matrix_quad_size);
    subtract_matrix(a11, a21, s9);

    matrix_t* s10 = create_matrix(matrix_quad_size);
    add_matrix(b11, b12, s10);

    matrix_t* p1 = create_matrix(matrix_quad_size);
    matrix_multiply_recursive_strassen(a11, s1, p1, k);

    matrix_t* p2 = create_matrix(matrix_quad_size);
    matrix_multiply_recursive_strassen(s2, b22, p2, k);

    matrix_t* p3 = create_matrix(matrix_quad_size);
    matrix_multiply_recursive_strassen(s3, b11, p3, k);

    matrix_t* p4 = create_matrix(matrix_quad_size);
    matrix_multiply_recursive_strassen(a22, s4, p4, k);

    matrix_t* p5 = create_matrix(matrix_quad_size);
    matrix_multiply_recursive_strassen(s5, s6, p5, k);

    matrix_t* p6 = create_matrix(matrix_quad_size);
    matrix_multiply_recursive_strassen(s7, s8, p6, k);

    matrix_t* p7 = create_matrix(matrix_quad_size);
    matrix_multiply_recursive_strassen(s9, s10, p7, k);

    add_matrix(p5, p4, c11);
    subtract_matrix(c11, p2, c11);
    add_matrix(c11, p6, c11);

    add_matrix(p1, p2, c12);

    add_matrix(p3, p4, c21);

    add_matrix(p5, p1, c22);
    subtract_matrix(c22, p3, c22);
    subtract_matrix(c22, p7, c22);

    copy_matrix_quadrants_out(c11, c12, c21, c22, c);

    delete_matrix(s1);
    delete_matrix(s2);
    delete_matrix(s3);
    delete_matrix(s4);
    delete_matrix(s5);
    delete_matrix(s6);
    delete_matrix(s7);
    delete_matrix(s8);
    delete_matrix(s9);
    delete_matrix(s10);

    delete_matrix(p1);
    delete_matrix(p2);
    delete_matrix(p3);
    delete_matrix(p4);
    delete_matrix(p5);
    delete_matrix(p6);
    delete_matrix(p7);

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

static void copy_matrix_quadrants(
    matrix_t* src,
    matrix_t* src11,
    matrix_t* src12,
    matrix_t* src21,
    matrix_t* src22)
{
    int final_rows = src->num_cols_rows;
    int final_cols = src->num_cols_rows;
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

static void copy_matrix_quadrants_out(
    matrix_t* src11,
    matrix_t* src12,
    matrix_t* src21,
    matrix_t* src22,
    matrix_t* dest)
{
    int final_rows = dest->num_cols_rows;
    int final_cols = dest->num_cols_rows;
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