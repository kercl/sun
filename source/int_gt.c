#include <string.h>
#include <stdlib.h>
#include <stdio.h> // TODO: dbg

#include "int_gt.h"
#include "sort.h"

/* Increments the first entry i0 of arr
 * for that arr[i0] < upper_lim[i0]
 * enforcing the constraints:
 * lower_lim <= arr <= upper_lim
 * and arr[i-1] <= arr[i]
 * 
 * arguments:
 *   arr       -- lower_lim <= arr <= upper_lim
 *   length    -- lenght of the arrays
 *   lower_lim -- lower boundary for arr
 *   upper_lim -- upper boundary for arr
 * 
 * returns:
 *   0 -- if the initial array fulfils
 *        arr < upper_lim
 *   1 -- if initial array fulfils
 *        arr = upper_lim
 */
int _array_increment_by_limits(int *arr, size_t length, int *lower_lim, int *upper_lim) {
    for(int i = 0; i < length; ++i) {
        if(arr[i] < upper_lim[i]) {
            ++arr[i];
            if(i > 0) {
                for(int j = 0; j < i; ++j)
                    arr[j] = arr[i];
            }
            return 0;
        }
    }
    memcpy(arr, lower_lim, sizeof(int) * length);
    return 1;
}

void gt_min_int_pattern(int *pattern, int *toprow, size_t length) {
    size_t n = (length * (length + 1)) >> 1,
           offset = 0;

    for(size_t i = 0, j = 0; i < n; ++i, ++j) {
        pattern[i] = toprow[offset + j];

        if(j == length - offset - 1) {
            offset++;
            j = -1;
        }
    }
}

void gt_min_int_pattern_transposed(int *pattern, int *toprow, size_t length) {
    int n = (length * (length + 1)) >> 1;

    size_t offset = 0;
    for(size_t i = 0; i < length; ++i) {
        for(size_t j = 0; j < length - i; ++j) {
            pattern[offset + j] = toprow[length - i - 1];
        }
        offset += length - i;
    }
}

int* gt_allocate_min_int_pattern(int *toprow, size_t length) {
    int n = (length * (length + 1)) >> 1;
    int *pattern = malloc(n * sizeof(int));

    gt_min_int_pattern(pattern, toprow, length);

    return pattern;
}

void gt_multi_transpose(int *patterns, int num_patterns, size_t length) {
    size_t m = length - 1, idx_a, idx_b,
           n = (length * (length + 1)) >> 1;
    int tmp;

    for(int i = 0; i < (length >> 1); ++i) {
        for(int j = i; j < length - i; ++j) {
            idx_a = GT_RL_IDX(  i,   j, length);
            idx_b = GT_RL_IDX(m-j, m-i, length);

            for(int k = 0; k < num_patterns * n; k += n) {
                tmp = patterns[k + idx_a];
                patterns[k + idx_a] = patterns[k + idx_b];
                patterns[k + idx_b] = tmp;
            }
        }
    }
}

void gt_transpose(int *pattern, size_t length) {
    gt_multi_transpose(pattern, 1, length);
}

/* Increments a transposed Gelfant-Tsetlin pattern
 * arguments:
 *   length -- Length of top row
 * 
 * returns:
 *   1 -- if pattern was incremented
 *   0 -- if pattern is maximal
 */
int _gt_increment_transposed(int *pattern_tr, int *min_pattern_tr, size_t length) {
    int carry = 1;
    size_t m = length - 1, offset;

    for(int i = 0; i < m; ++i) {
        offset = GT_ROW_OFFSET(i, length);
        carry = _array_increment_by_limits(
            pattern_tr + offset, 
            m - i,
            min_pattern_tr + offset,
            pattern_tr + offset + length - i);
        
        if(carry) {          
            if(i < m - 1)  
                memcpy(pattern_tr, min_pattern_tr, offset + length - i);
            else
                return 0;
        } else
            return 1;
    }
    return 0;
}

size_t _ggt_pos(size_t x, size_t y) {
  size_t c;
  while (y != 0) {
    c = x % y;
    x = y;
    y = c;
  }
  return x;
}

size_t gt_num_of_patterns(int *toprow, size_t length) {
    int *m = malloc(sizeof(int) * (length - 1));

    for(int i = 0; i < length - 1; ++i)
        m[i] = toprow[i] - toprow[i + 1];

    size_t num = 1, denom = 1,
           red_factor, factor;
    for(int i = 1; i <= length - 1; ++i) {
        for(int j = 0; j < length - i; ++j) {
            factor = 0;

            for(int k = 0; k < i; ++k)
                factor += m[j + k];
            num *= factor + i;
            denom *= i;

            red_factor = _ggt_pos(num, denom);
            num /= red_factor;
            denom /= red_factor;
        }
    }

    return num;
}

void gt_generate_all_transposed(int **patterns, size_t *num_patterns, int *toprow, size_t length) {
    size_t n_patterns = gt_num_of_patterns(toprow, length),
           n_entries = (length * (length + 1)) >> 1;
    int *_patterns = malloc(sizeof(int) * n_entries * n_patterns);

    gt_min_int_pattern_transposed(_patterns, toprow, length);

    for(int i = 1; i < n_patterns; ++i) {
        memcpy(_patterns + i * n_entries,
               _patterns + (i - 1) * n_entries,
               n_entries * sizeof(int));
        _gt_increment_transposed(_patterns + i * n_entries, _patterns, length);
    }

    *patterns = _patterns;
    *num_patterns = n_patterns;
}

void gt_generate_all(int **patterns, size_t *num_patterns, int *toprow, size_t length) {
    gt_generate_all_transposed(patterns, num_patterns, toprow, length);
    gt_multi_transpose(*patterns, *num_patterns, length);
}

void gt_sort_patterns(int *patterns, size_t n_patterns, size_t length) {
    sort_pattern(patterns, n_patterns, length);
}
