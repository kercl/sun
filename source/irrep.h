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

#ifndef SOURCE_IRREP_H_
#define SOURCE_IRREP_H_

#include <stdint.h>

#include "int_gt.h"

typedef int16_t mat_int_t;

/* Construct the l-th center generator
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
 */
void gt_center_generator_diag(struct gt_tree *patterns,
                              int l,
                              mat_int_t *diagonal);

/* Construct the l-th lowering operator
 * For GT-patterns of length n, note that
 *
 * l = 1,...,n-1
 * 
 * TODO: add documentation
 */
void gt_lowering_operator(struct gt_tree *patterns, int l);

#endif  // SOURCE_IRREP_H_
