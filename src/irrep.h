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

#ifndef SUN_CORE_IRREP_H_
#define SUN_CORE_IRREP_H_

#include <stdint.h>

#include "int_gt.h"

typedef int mat_int_t;

/* Convert dynkin label to GT top row
 */
gt_int_t*
gt_top_row_from_dynkin(gt_int_t *dynkin, size_t length);

/* Construct the l-th generator of the Cartan
 * subalgebra.
 * For GT-patterns of length n, note that
 *
 * l = 1,...,n-1
 * 
 * numerator and denominator must be arrays
 * of length tree->num_patterns. The entries
 * are the diagonal entries of the center
 * generators.
 * Note: the matrices are normed s.t. the
 * diagonals are integer valued.
 * To obtain the same entries as in [ref]
 * divide by 2.
 */
void
csa_generator_diag_from_gt(struct gt_tree *patterns,
                           size_t l,
                           mat_int_t *diagonal);

/* Construct the l-th lowering operator
 * For GT-patterns of length n, note that
 *
 * l = 1,...,n-1
 * 
 * TODO: add documentation
 */
size_t
lowering_operator_from_gt(struct gt_tree *patterns,
                          size_t l,
                          mat_int_t **ptr_numerators,
                          mat_int_t **ptr_denominators,
                          size_t **ptr_row,
                          size_t **ptr_col);

/* Compute dimension from a given top row
 */
size_t
dimension_from_dynkin(gt_int_t *dynkin, size_t length);

#endif  // SUN_CORE_IRREP_H_
