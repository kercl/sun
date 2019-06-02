#include "irrep.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void gt_center_generator_diag(struct gt_tree *patterns,
                              int l,
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

void gt_lowering_operator(struct gt_tree *patterns, int l) {
    size_t length = patterns->length,
           row_len = patterns->length - l,
           num_patterns = patterns->num_patterns;
    gt_int_t *pattern_array = patterns->array_representation;

    size_t n_entries = (length * (length + 1)) / 2,
           row_start = n_entries - (row_len * (row_len + 1)) / 2,
           Mj;

    mat_int_t *numerators = malloc(sizeof(mat_int_t) * num_patterns),
              *denominators = malloc(sizeof(mat_int_t) * num_patterns);
    int8_t *imaginary = malloc(sizeof(int8_t) * num_patterns);

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
