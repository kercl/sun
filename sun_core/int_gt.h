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

#ifndef SUN_CORE_INT_GT_H_
#define SUN_CORE_INT_GT_H_

#include <stddef.h>
#include <stdint.h>

typedef int16_t gt_int_t;

/*
 * Simplification of r*n+((r-1)*r)/2+c-r
 */
#define GT_RL_IDX(r, c, n) ((((c) << 1) + (((n) << 1) - (r) - 1) * (r)) >> 1)
#define GT_ROW_OFFSET(r, n) ((r) * (n) - ((((r) - 1) * (r)) >> 1))

struct gt_node {
    struct gt_node *children;
    size_t index;
};

struct gt_tree {
    gt_int_t *array_representation;
    struct gt_node root;
    size_t length, num_patterns;
};

/* Transpose a given pattern of length `length`.
 * A transposed pattern is defined by
 * 
 * m11 ... m1N           mN1 ... m1N
 *  .       .             .       .
 *   .     .      --->     .     .
 *    .   .                 .   .
 *     mN1                   m11
 * 
 * Example:
 * 
 * 2 1 0        0 0 0
 *  1 0   --->   1 1
 *   0            2
 */
void
gt_transpose(gt_int_t *pattern, size_t length);

/* Compute the transposed of multiple patterns
 * For the definition of a transposed pattern
 * see `gt_transpose`
 */
void
gt_multi_transpose(gt_int_t *patterns, size_t num_patterns, size_t length);

/* For a given top row, allocate a 
 * new pattern and fill it with the unique
 * entries such that the sum of those
 * entries is minimal among all patterns
 * sharing the same top row.
 * 
 * Example: toprow=(2,1,0)
 * 
 * Minimal pattern:
 * 
 * 2 1 0
 *  1 0
 *   0
 * 
 * Result:
 * (2,1,0,1,0,0)
 * 
 * arguments:
 *   toprow -- array containing the top row
 *             constraint
 *   length -- length of top row
 * 
 * returns:
 *   pointer to pattern array
 */
gt_int_t*
gt_allocate_min_int_pattern(gt_int_t *toprow, size_t length);

/* Generate all possible GT patterns for
 * a given top-row.
 */
gt_int_t*
gt_patterns_generate(gt_int_t *toprow, size_t length);

/* Generate all possible GT patterns for
 * a given top row in transposed form
 */
void
gt_generate_all_transposed(gt_int_t **pattern,
                           size_t *num_entries,
                           gt_int_t *toprow,
                           size_t length);

/* Generate all possible GT patterns for
 * a given top row in transposed form
 */
void
gt_generate_all(gt_int_t **pattern,
                size_t *num_entries,
                gt_int_t *toprow,
                size_t length);

/* Calculate number of all possible patterns
 * for a given top row by means of Weyls formula
 */
size_t
gt_num_of_patterns(gt_int_t *toprow, size_t length);

/* Generates a search tree from a list of
 * patterns and store the indices of that list
 * at the leaf nodes.
 * 
 * Example: Consider the list of patterns
 * 
 * 2 0 0  2 0 0  2 0 0  2 0 0  2 0 0  2 0 0
 *  0 0    1 0    1 0    2 0    2 0    2 0
 *   0      0      1      0      1      2
 * 
 *      .--(0)--[ ]--(0)--[ ]---(0)--[idx: 0]
 *     /
 *    /                      .--(0)--[idx: 1]
 *   /                      /
 * [ ]-----(1)--[ ]--(0)--[ ]---(1)--[idx: 2]
 *   \
 *    \                      .--(0)--[idx: 3]
 *     \                    /
 *      '--(2)--[ ]--(0)--[ ]---(1)--[idx: 4]
 *                          \
 *                           '--(2)--[idx: 5]
 * 
 * Storing the patterns in this structure helps
 * improving the retrieval of the matrix indices
 * for the raising operators.
 * 
 * Note that `array_representation` is set to
 * `patterns` and is used int the computation
 * of the matrix representations. It should 
 * be deallocated only after `gt_free_tree`.
 */
void
gt_list_to_tree(struct gt_tree *tree,
                gt_int_t *patterns,
                size_t num_patterns,
                size_t length);

/* Free tree structure generated by gt_list_to_tree
 * Note that it only deallocates the nodes in the 
 * tree and not the list of patterns `array_representation`.
 */
void
gt_free_tree(struct gt_tree *tree);

/* For a given pattern, find the corresponding indes
 * in `tree->array_representation`.
 */
size_t
gt_locate_in_tree(struct gt_tree *tree, gt_int_t *pattern);

#endif  // SUN_CORE_INT_GT_H_
