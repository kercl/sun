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

#include "sun_core/irrep.h"

#include <stdlib.h>
#include <stdint.h>

gt_int_t*
gt_top_row_from_dynkin(gt_int_t *dynkin, size_t length) {
    gt_int_t *top_row = malloc(sizeof(gt_int_t) * (length + 1));

    top_row[length] = 0;
    for (int i = length - 1; i >= 0; i--)
        top_row[i] = top_row[i + 1] + dynkin[i];

    return top_row;
}

void
center_generator_diag_from_gt(struct gt_tree *patterns,
                                   size_t l,
                                   mat_int_t *diagonal) {
    size_t length = patterns->length,
           row_len = patterns->length - l,
           num_patterns = patterns->num_patterns;
    gt_int_t *pattern_array = patterns->array_representation;

    size_t n_entries = (length * (length + 1)) / 2,
           row_start = n_entries - (row_len * (row_len + 1)) / 2;

    for (size_t M = 0, Mi = 0; Mi < num_patterns; M+=n_entries, Mi++) {
        diagonal[Mi] = 0;
        for (size_t j = 0; j < row_len; j++) {
            if (j < row_len) {
                // sum over row l
                diagonal[Mi] +=
                    2 * pattern_array[M + row_start + j];
            }
            if (j < row_len - 1 && l < length) {
                // sum over row l-1
                diagonal[Mi] -=
                    pattern_array[M + row_start + row_len + j];
            }
            // sum over row l+1
            diagonal[Mi] -=
                pattern_array[M + row_start - (row_len + 1) + j];
        }
    }
}

void
lowering_operator_from_gt(struct gt_tree *patterns, size_t l) {
    size_t length = patterns->length,
           row_len = patterns->length - l,
           num_patterns = patterns->num_patterns;
    gt_int_t *pattern_array = patterns->array_representation;

    size_t n_entries = (length * (length + 1)) / 2,
           row_start = n_entries - (row_len * (row_len + 1)) / 2;
    int Mj;

    // mat_int_t *numerators = malloc(sizeof(mat_int_t) * num_patterns),
    //           *denominators = malloc(sizeof(mat_int_t) * num_patterns);
    // int8_t *imaginary = malloc(sizeof(int8_t) * num_patterns);

    int M_lk, M_lkp;
    for (size_t M = 0, Mi = 0; Mi < num_patterns; M+=n_entries, Mi++) {
        for (size_t k = 0; k < row_len; k++) {
            // decrement pattern M -> M - M^{k,l}
            pattern_array[M + row_start + k]--;
            Mj = gt_locate_in_tree(patterns, pattern_array + M);
            // recover pattern
            pattern_array[M + row_start + k]++;

            // pattern does not exist
            if (Mj == -1)
                continue;

            int numerator = 1, denominator = 1;

            for (size_t k_p = 0; k_p < row_len + 1 && numerator != 0; k_p++) {
                M_lk = pattern_array[M + row_start + k];

                if (k_p < row_len && k != k_p) {
                    M_lkp = pattern_array[M + row_start + k_p];
                    denominator *=
                        (M_lkp - M_lk + k - k_p + 1) *
                        (M_lkp - M_lk + k - k_p);
                }
                if (k_p < row_len - 1) {
                    numerator *= pattern_array[M + row_start + row_len + k_p]
                                 - M_lk + k - k_p;
                }
                numerator *= pattern_array[M + row_start - (row_len + 1) + k_p]
                             - M_lk + k - k_p + 1;
            }
        }
    }
}

size_t
dimension_from_dynkin(gt_int_t *dynkin, size_t length) {
    gt_int_t *top_row = gt_top_row_from_dynkin(dynkin, length);
    size_t dim = gt_num_of_patterns(top_row, length + 1);
    free(top_row);
    return dim;
}
