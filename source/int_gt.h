/* Library for Gelfant-Tsetlin pattern generation
 * for SU(n) representation generator.
 */

#ifndef GT_H
#define GT_H

#include <stddef.h>

/*
 * Simplification of r*n+((r-1)*r)/2+c-r
 */
#define GT_RL_IDX(r,c,n) ((((c) << 1) + (((n) << 1) - (r) - 1) * (r)) >> 1)
#define GT_ROW_OFFSET(r,n) ((r) * (n) - ((((r) - 1) * (r)) >> 1))

/*
 * TODO: documentation
 */
void gt_transpose(int *pattern, size_t length);

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
int* gt_allocate_min_int_pattern(int *toprow, size_t length);

/* Generate all possible GT patterns for
 * a given top-row.
 */
int* gt_patterns_generate(int *toprow, size_t length);

/* Generate all possible GT patterns for
 * a given top row in transposed form
 */
void gt_generate_all_transposed(int **pattern, int *num_entries, int *toprow, size_t length);

/* Generate all possible GT patterns for
 * a given top row in transposed form
 */
void gt_generate_all(int **pattern, int *num_entries, int *toprow, size_t length);

/* Calculate number of all possible patterns
 * for a given top row.
 */
size_t gt_num_of_patterns(int *toprow, size_t length);

/* Sort patterns lexicographically for quick
 * lookups.
 */
void gt_sort_patterns(int *patterns, size_t num_entries, size_t length);

/*
 *
 */
size_t gt_index_of(int *patterns, size_t num_entries, size_t length, size_t idx);

#endif
