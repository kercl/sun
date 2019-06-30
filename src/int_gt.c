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

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "int_gt.h"

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
uint8_t
_array_increment_by_limits(gt_int_t *arr,
                           size_t length,
                           gt_int_t *lower_lim,
                           gt_int_t *upper_lim) {
    for (size_t i = 0; i < length; ++i) {
        if (arr[i] < upper_lim[i]) {
            ++arr[i];
            if (i > 0) {
                for (size_t j = 0; j < i; ++j)
                    arr[j] = arr[i];
            }
            return 0;
        }
    }
    memcpy(arr, lower_lim, sizeof(gt_int_t) * length);
    return 1;
}

void
gt_min_int_pattern(gt_int_t *pattern, gt_int_t *toprow, size_t length) {
    size_t n = (length * (length + 1)) >> 1,
           offset = 0;

    for (size_t i = 0, j = 0; i < n; i++) {
        pattern[i] = toprow[offset + j];

        if (j == length - offset - 1) {
            offset++;
            j = 0;
        } else {
            j++;
        }
    }
}

void
gt_min_int_pattern_transposed(gt_int_t *pattern,
                              gt_int_t *toprow,
                              size_t length) {
    size_t offset = 0;
    for (size_t i = 0; i < length; ++i) {
        for (size_t j = 0; j < length - i; ++j) {
            pattern[offset + j] = toprow[length - i - 1];
        }
        offset += length - i;
    }
}

gt_int_t*
gt_allocate_min_int_pattern(gt_int_t *toprow, size_t length) {
    size_t n = (length * (length + 1)) >> 1;
    gt_int_t *pattern = malloc(n * sizeof(gt_int_t));

    gt_min_int_pattern(pattern, toprow, length);

    return pattern;
}

void
gt_multi_transpose(gt_int_t *patterns, size_t num_patterns, size_t length) {
    size_t m = length - 1, idx_a, idx_b,
           n = (length * (length + 1)) >> 1;
    int tmp;

    for (size_t i = 0; i < (length >> 1); ++i) {
        for (size_t j = i; j < length - i; ++j) {
            idx_a = GT_RL_IDX(i, j, length);
            idx_b = GT_RL_IDX(m-j, m-i, length);

            for (size_t k = 0; k < num_patterns * n; k += n) {
                tmp = patterns[k + idx_a];
                patterns[k + idx_a] = patterns[k + idx_b];
                patterns[k + idx_b] = tmp;
            }
        }
    }
}

void
gt_transpose(gt_int_t *pattern, size_t length) {
    gt_multi_transpose(pattern, 1, length);
}

/* Increments a transposed Gelfand-Tsetlin pattern
 * arguments:
 *   length -- Length of top row
 * 
 * returns:
 *   1 -- if pattern was incremented
 *   0 -- if pattern is maximal
 */
int
_gt_increment_transposed(gt_int_t *pattern_tr,
                         gt_int_t *min_pattern_tr,
                         size_t length) {
    uint8_t carry = 1;
    size_t m = length - 1, offset;

    for (size_t i = 0; i < m; ++i) {
        offset = GT_ROW_OFFSET(i, length);
        carry = _array_increment_by_limits(
            pattern_tr + offset,
            m - i,
            min_pattern_tr + offset,
            pattern_tr + offset + length - i);

        if (carry) {
            if (i < m - 1)
                memcpy(pattern_tr, 
                       min_pattern_tr, 
                       sizeof(gt_int_t) * (offset + length - i));
            else
                return 0;
        } else {
            return 1;
        }
    }
    return 0;
}

size_t
_abs_gcd(size_t x, size_t y) {
  size_t c;
  while (y != 0) {
    c = x % y;
    x = y;
    y = c;
  }
  return x;
}

size_t
gt_num_of_patterns(gt_int_t *toprow, size_t length) {
    gt_int_t *m = malloc(sizeof(gt_int_t) * (length - 1));

    for (size_t i = 0; i < length - 1; ++i)
        m[i] = toprow[i] - toprow[i + 1];

    size_t num = 1, denom = 1,
           red_factor, factor;
    for (size_t i = 1; i <= length - 1; ++i) {
        for (size_t j = 0; j < length - i; ++j) {
            factor = 0;

            for (size_t k = 0; k < i; ++k)
                factor += m[j + k];
            num *= factor + i;
            denom *= i;

            red_factor = _abs_gcd(num, denom);
            num /= red_factor;
            denom /= red_factor;
        }
    }

    return num;
}

void
gt_generate_all_transposed(gt_int_t **patterns,
                           size_t *num_patterns,
                           gt_int_t *toprow,
                           size_t length) {
    size_t n_patterns = gt_num_of_patterns(toprow, length),
           n_entries = (length * (length + 1)) >> 1;
    gt_int_t *_patterns = malloc(sizeof(gt_int_t) * n_entries * n_patterns);

    gt_min_int_pattern_transposed(_patterns, toprow, length);

    for (size_t i = 1; i < n_patterns; ++i) {
        memcpy(_patterns + i * n_entries,
               _patterns + (i - 1) * n_entries,
               n_entries * sizeof(gt_int_t));
        _gt_increment_transposed(_patterns + i * n_entries, _patterns, length);
    }

    *patterns = _patterns;
    *num_patterns = n_patterns;
}

void
gt_generate_all(gt_int_t **patterns,
                size_t *num_patterns,
                gt_int_t *toprow,
                size_t length) {
    gt_generate_all_transposed(patterns, num_patterns, toprow, length);
    gt_multi_transpose(*patterns, *num_patterns, length);
}

void
_insert_pattern(struct gt_node *node,
                gt_int_t *pattern,
                size_t index,
                size_t length) {
    /* patterns always start with the same length
    entries. We don't have to store them to the
    tree. */

    // we only allow patterns that contain more than just the top row
    assert(length >= 2);

    gt_int_t *toprow = pattern;
    size_t d_limit = length - 2;

    for (size_t d = length, i = 0; d < (length * (length + 1)) >> 1; d++) {
        if (node->children == NULL) {
            node->children = malloc(sizeof(struct gt_node) * (toprow[i] + 1));
            for (int j = 0; j <= toprow[i]; ++j) {
                node->children[j].children = NULL;
                node->children[j].index = -1;
            }
        }
        assert(pattern[d] <= pattern[i]);
        if (i == d_limit) {
            i = 0;
            d_limit--;
        } else {
            i++;
        }

        node = node->children + pattern[d] - toprow[length - 1];
    }
    node->index = index;
}

size_t
gt_locate_in_tree(struct gt_tree *tree, gt_int_t *pattern) {
    size_t length = tree->length;

    // we only allow patterns that contain more than just the top row
    assert(length >= 2);

    gt_int_t *toprow = pattern;
    size_t d_limit = length - 2;

    struct gt_node *node = &tree->root;

    for (size_t d = length, i = 0; d < (length * (length + 1)) >> 1; d++) {
        if (node->children == NULL
            || pattern[d] > pattern[i]
            || pattern[d] < toprow[length - 1])
            return -1;

        if (i == d_limit) {
            i = 0;
            d_limit--;
        } else {
            i++;
        }

        node = node->children + pattern[d] - toprow[length - 1];
    }
    return node->index;
}

void
gt_list_to_tree(struct gt_tree *tree,
                gt_int_t *patterns,
                size_t num_patterns,
                size_t length) {
    tree->array_representation = patterns;
    tree->length = length;
    tree->num_patterns = num_patterns;
    tree->root.children = NULL;

    size_t m = (length * (length + 1)) >> 1;

    for (size_t i = 0; i < num_patterns; i++) {
        _insert_pattern(&tree->root, patterns + i * m, i, length);
    }
}

void
_free_node(gt_int_t *toprow, size_t d, size_t length, struct gt_node *node) {
    if (length == 0 || node->children == NULL)
        return;

    if (d == length) {
        d = 0;
        length--;
    }

    for (gt_int_t i = 0; i <= toprow[d]; i++)
        _free_node(toprow, d+1, length, node->children + i);

    free(node->children);
}

void
gt_free_tree(struct gt_tree *tree, int free_array) {
    _free_node(tree->array_representation, 0, tree->length - 1, &tree->root);
    if(free_array)
        free(tree->array_representation);
    tree->num_patterns = 0;
}
