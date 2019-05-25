/* Library for Gelfant-Tsetlin pattern generation
 * for SU(n) representation generator.
 */

#ifndef GT_H
#define GT_H

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
int _array_increment_by_limits(int *arr, size_t length, int *lower_lim, int *upper_lim);

/* Increments a Gelfant-Tsetlin pattern
 * arguments:
 *   length -- Length of top row
 * 
 * returns:
 *   1 -- if pattern was incremented
 *   0 -- if pattern is maximal
 */
int _gt_increment(int *pattern, int *min_pattern, size_t length);

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
int* gt_min_int_pattern(int *toprow, size_t length);

#endif
