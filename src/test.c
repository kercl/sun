/* Copyright (c) 2018 Clemens Kerschbaumer
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "int_gt.h"
#include "irrep.h"

int _array_increment_by_limits(gt_int_t *arr,
                               size_t length,
                               gt_int_t *lower_lim,
                               gt_int_t *upper_lim);

int _gt_increment_transposed(gt_int_t *pattern_tr,
                             gt_int_t *min_pattern_tr,
                             size_t length);

void print_pattern_raligned(gt_int_t *pattern, size_t length) {
    for (int i = 0; i < length; ++i) {
        for (int k = 0; k < i; k++)
            printf("   ");
        for (int j = i; j < length; ++j) {
            printf("%2d ", pattern[GT_RL_IDX(i, j, length)]);
        }
        printf("\n");
    }
}

void print_pattern_flattened(int *pattern, size_t length) {
    for (int i = 0; i < (length*(length + 1)) >> 1; ++i) {
        printf("%2d ", pattern[i]);
    }
}

void test__array_increment_by_limits() {
    gt_int_t test_data[48][5] = {
     {4, 4, 1, 1, 0}, {2, 2, 2, 1, 0}, {3, 2, 2, 1, 0},
     {4, 2, 2, 1, 0}, {3, 3, 2, 1, 0}, {4, 3, 2, 1, 0},
     {4, 4, 2, 1, 0}, {3, 3, 3, 1, 0}, {4, 3, 3, 1, 0},
     {4, 4, 3, 1, 0}, {2, 2, 2, 2, 0}, {3, 2, 2, 2, 0},
     {4, 2, 2, 2, 0}, {3, 3, 2, 2, 0}, {4, 3, 2, 2, 0},
     {4, 4, 2, 2, 0}, {3, 3, 3, 2, 0}, {4, 3, 3, 2, 0},
     {4, 4, 3, 2, 0}, {1, 1, 1, 1, 1}, {2, 1, 1, 1, 1},
     {3, 1, 1, 1, 1}, {4, 1, 1, 1, 1}, {2, 2, 1, 1, 1},
     {3, 2, 1, 1, 1}, {4, 2, 1, 1, 1}, {3, 3, 1, 1, 1},
     {4, 3, 1, 1, 1}, {4, 4, 1, 1, 1}, {2, 2, 2, 1, 1},
     {3, 2, 2, 1, 1}, {4, 2, 2, 1, 1}, {3, 3, 2, 1, 1},
     {4, 3, 2, 1, 1}, {4, 4, 2, 1, 1}, {3, 3, 3, 1, 1},
     {4, 3, 3, 1, 1}, {4, 4, 3, 1, 1}, {2, 2, 2, 2, 1},
     {3, 2, 2, 2, 1}, {4, 2, 2, 2, 1}, {3, 3, 2, 2, 1},
     {4, 3, 2, 2, 1}, {4, 4, 2, 2, 1}, {3, 3, 3, 2, 1},
     {4, 3, 3, 2, 1}, {4, 4, 3, 2, 1}, {3, 3, 1, 1, 0}};

    gt_int_t arr[] = {4, 3, 1, 1, 0},
             lower[] = {3, 3, 1, 1, 0},
             upper[] = {4, 4, 3, 2, 1};

    int res = 0;
    for (int i = 0; i < 48; i++) {
        assert(res == 0);
        res = _array_increment_by_limits(arr, 5, lower, upper);
        assert(memcmp(arr, test_data[i], sizeof(gt_int_t) * 5) == 0);
    }
    assert(res == 1);
}

void test_gt_allocate_min_int_pattern() {
    gt_int_t toprow[]  = {5, 5, 3, 1, 0},
             toprow2[] = {9, 7, 1, 0};

    gt_int_t result[]  = {5, 5, 3, 1, 0, 5, 3, 1, 0, 3, 1, 0, 1, 0, 0},
             result2[] = {9, 7, 1, 0, 7, 1, 0, 1, 0, 0};

    gt_int_t *pattern = gt_allocate_min_int_pattern(toprow, 5);
    assert(memcmp(pattern, result, 15 * sizeof(gt_int_t)) == 0);
    free(pattern);

    pattern = gt_allocate_min_int_pattern(toprow2, 4);
    assert(memcmp(pattern, result2, 10 * sizeof(gt_int_t)) == 0);
    free(pattern);
}

void test_gt_transpose() {
    gt_int_t tr_1[] = {6, 5, 4, 3, 2, 1, 0},
             tr_2[] = {8, 7, 6, 5, 4};

    gt_int_t res_1[] = {
         0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1,
               2, 2, 2, 2, 2,
                  3, 3, 3, 3,
                     4, 4, 4,
                        5, 5,
                           6};
    gt_int_t res_2[] = {
         4, 4, 4, 4, 4,
            5, 5, 5, 5,
               6, 6, 6,
                  7, 7,
                     8};

    gt_int_t *pattern = gt_allocate_min_int_pattern(tr_1, 7);
    gt_transpose(pattern, 7);
    assert(memcmp(pattern, res_1, 28*sizeof(gt_int_t)) == 0);
    free(pattern);

    pattern = gt_allocate_min_int_pattern(tr_2, 5);
    gt_transpose(pattern, 5);
    assert(memcmp(pattern, res_2, 15*sizeof(gt_int_t)) == 0);
    free(pattern);
}

void test_gt_generate_all_transposed_1() {
    gt_int_t toprow[] = {2, 1, 0};
    size_t length = 3;

    gt_int_t *patterns;
    size_t n_entries;

    assert(gt_num_of_patterns(toprow, length) == 8);

    gt_generate_all(&patterns, &n_entries, toprow, length);
    struct gt_tree tree;
    gt_list_to_tree(&tree, patterns, n_entries, length);

    int m = length * (length + 1) >> 1;

    for (int i = 0; i < n_entries; ++i)
        assert(i == gt_locate_in_tree(&tree, patterns + i * m));
    assert(n_entries == 8);

    gt_free_tree(&tree, 1);
}

void test_gt_generate_all_transposed_2() {
    gt_int_t toprow[] = {3, 2, 1, 0};
    size_t length = 4;

    gt_int_t *patterns;
    size_t n_entries;

    assert(gt_num_of_patterns(toprow, length) == 64);

    gt_generate_all(&patterns, &n_entries, toprow, length);
    struct gt_tree tree;
    gt_list_to_tree(&tree, patterns, n_entries, length);

    int m = length * (length + 1) >> 1;

    for (int i = 0; i < n_entries; ++i)
        assert(i == gt_locate_in_tree(&tree, patterns + i * m));
    assert(n_entries == 64);

    gt_free_tree(&tree, 1);
}

void test_dimension_from_dynkin() {
    gt_int_t dynkin_1[] = {1},
             dynkin_2[] = {1, 1},
             dynkin_3[] = {1, 0},
             dynkin_4[] = {1, 0, 0};

    assert(dimension_from_dynkin(dynkin_1, 1) == 2);
    assert(dimension_from_dynkin(dynkin_2, 2) == 8);
    assert(dimension_from_dynkin(dynkin_3, 2) == 3);
    assert(dimension_from_dynkin(dynkin_4, 3) == 4);
}

int main(int argc, char **argv) {
    test__array_increment_by_limits();
    test_gt_allocate_min_int_pattern();
    test_gt_transpose();
    test_gt_generate_all_transposed_2();
    test_dimension_from_dynkin();

    /*gt_int_t toprow[] = {2, 1, 0};
    size_t length = sizeof(toprow) / sizeof(gt_int_t);

    gt_int_t *patterns;
    size_t n_patterns;

    gt_generate_all(&patterns, &n_patterns, toprow, length);
    struct gt_tree tree;
    gt_list_to_tree(&tree, patterns, n_patterns, length);

    for (int i = 0; i < n_patterns; i++)
        print_pattern_raligned(patterns
            + i * (length * (length + 1) / 2), length);

    mat_int_t *x3 = malloc(sizeof(mat_int_t) * n_patterns);

    gt_center_generator_diag(&tree, 1, x3);

    for (int i = 0; i < n_patterns; i++) {
        printf("%d/2 ", x3[i]);
    }
    printf("\n");

    printf("--\n");
    gt_lowering_operator(&tree, 2);
    printf("--\n");
    gt_lowering_operator(&tree, 1);

    printf("\n");*/
}
