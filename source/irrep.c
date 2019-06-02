#include "irrep.h"
#include "sort.h"

#include <stdio.h>

void _sigma(int *result, int *patterns, size_t n_patterns, size_t length, int l) {
    size_t n_entries = (length * (length + 1)) / 2;
    size_t row_start = n_entries - ((length - l) * (length - l + 1)) / 2;
    
    for(size_t i = 0; i < n_patterns * n_entries; i+=n_entries) {
        result[i] = 0;
        for(size_t j = 0; j < length - l; j++)
            result[i] += patterns[i + row_start + j];
    }
}

/*
 * l = 1,...,length - 1
 */
void gt_center_generator(struct gt_tree *patterns, int l, int *numerator, int *denominator) {
    size_t length = patterns->length,
           row_length = patterns->length - l,
           num_patterns = patterns->num_patterns;
    int *pattern_array = patterns->array_representation;

    size_t n_entries = (length * (length + 1)) / 2,
           row_start = n_entries - (row_length * (row_length + 1)) / 2;

    for(size_t M = 0, Mi = 0; Mi < num_patterns; M+=n_entries, Mi++) {
        numerator[Mi] = 0;
        for(size_t j = 0; j < row_length; j++) {
            if(j < row_length) {
                numerator[Mi] += 2 * pattern_array[M + row_start + j]; // sum over row l
            }
            if(j < row_length - 1 && l < length) {
                numerator[Mi] -= pattern_array[M + row_start + row_length + j]; // sum over row l-1
            }
            numerator[Mi] -= pattern_array[M + row_start - (row_length + 1) + j]; // sum over row l+1
        }
        denominator[Mi] = 2;
    }
}

/*
 * l = 1,...,length - 1
 */
void gt_lowering_operator(struct gt_tree *patterns, int l) {
    size_t length = patterns->length,
           row_length = patterns->length - l,
           num_patterns = patterns->num_patterns;
    int *pattern_array = patterns->array_representation;

    size_t n_entries = (length * (length + 1)) / 2,
           row_start = n_entries - (row_length * (row_length + 1)) / 2,
           Mj;

    for(size_t M = 0, Mi = 0; Mi < num_patterns; M+=n_entries, Mi++) {
        for(size_t k = 0; k < row_length; k++) {
            pattern_array[M + row_start + k]--; // decrement pattern M -> M - M^{k,l}
            Mj = gt_locate_in_tree(patterns, pattern_array + M);
            pattern_array[M + row_start + k]++; // recover pattern

            if(Mj == -1) // pattern does not exist
                continue;

            int numerator = 1, denominator = 1;

            for(size_t k_prime = 0; k_prime < row_length + 1; k_prime++) {
                int M_lk = pattern_array[M + row_start + k];

                if(k_prime < row_length && k != k_prime) {
                    denominator *= (pattern_array[M + row_start + k_prime] - M_lk + k - k_prime + 1) * 
                                   (pattern_array[M + row_start + k_prime] - M_lk + k - k_prime);
                }
                if(k_prime < row_length - 1) {
                    printf(" * (%d)", pattern_array[M + row_start + row_length + k_prime] - M_lk + k - k_prime);
                    numerator *= pattern_array[M + row_start + row_length + k_prime] - M_lk + k - k_prime;
                }
                //if(l < length - 1) {
                printf(" * %d", pattern_array[M + row_start - (row_length + 1) + k_prime] - M_lk + k - k_prime + 1);
                numerator *= pattern_array[M + row_start - (row_length + 1) + k_prime] - M_lk + k - k_prime + 1;
                //}
            }
            printf("\n");

            printf("Matrix element %ld,%ld --> %d / %d\n", Mj, Mi, numerator, denominator);
        }
    }
}
