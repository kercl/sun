#include "irrep.h"

void _sigma(int *result, int *patterns, size_t n_patterns, size_t length, int l) {
    size_t n_entries = (length * (length + 1)) / 2;
    size_t row_start = n_entries - ((length - l) * (length - l + 1)) / 2;
    
    for(size_t i = 0; i < n_patterns * n_entries; i+=n_entries) {
        result[i] = 0;
        for(size_t j = 0; j < length - l; j++)
            result[i] += patterns[i + row_start + j];
    }
}

void _add(int *result, int *a, int *b, size_t length) {
    for(size_t i = 0; i < length; i++) {
        result[i] = a[i] + b[i];
    }
}

void gt_center_generator(struct gt_tree *patterns, int l, int *numerator, int *denominator) {
    size_t length = patterns->length,
           num_patterns = patterns->num_patterns;
    int *pattern_array = patterns->array_representation;

    size_t n_entries = (length * (length + 1)) / 2;
    size_t row_start = n_entries - ((length - l) * (length - l + 1)) / 2;

    for(size_t M = 0, Mi; M < num_patterns * n_entries; M+=n_entries, Mi++) {
        numerator[Mi] = 0;
        for(size_t j = 0; j < length - l + 1; j++) {
            if(j < length - l) {
                numerator[Mi] += 2 * pattern_array[M + row_start + j];
            }
            if(j < length - l - 1 && l < length) {
                numerator[Mi] -= pattern_array[M + row_start + j + length - l];
            }
            numerator[Mi] -= pattern_array[M + row_start + j - (length - l) - 1];
        }
        denominator[Mi] = 2;
    }
}
