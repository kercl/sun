#include <string.h>
#include <stdlib.h>
#include <stdio.h> # TODO: remove, dbg

#include "int_gt.h"

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

int* gt_min_int_pattern(int *toprow, size_t length) {
    int n = length * (length + 1) >> 1;
    int *pattern = malloc(n * sizeof(int));

    int oinc = length;

    for(int i = 0, j = 0; i < n; ++i, ++j) {
        pattern[i] = toprow[j];

        oinc--;
        if(oinc == 0) {
            oinc--;
            j = -1;
        }
    }
    return pattern;
}

int _gt_increment(int *pattern, int *min_pattern, size_t length) {
    
}