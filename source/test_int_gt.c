#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "int_gt.h"

void test__array_increment_by_limits() {
    int test_data[48][5] = {
     {4,4,1,1,0}, {2,2,2,1,0}, {3,2,2,1,0}, {4,2,2,1,0}, {3,3,2,1,0}, {4,3,2,1,0},
     {4,4,2,1,0}, {3,3,3,1,0}, {4,3,3,1,0}, {4,4,3,1,0}, {2,2,2,2,0}, {3,2,2,2,0},
     {4,2,2,2,0}, {3,3,2,2,0}, {4,3,2,2,0}, {4,4,2,2,0}, {3,3,3,2,0}, {4,3,3,2,0},
     {4,4,3,2,0}, {1,1,1,1,1}, {2,1,1,1,1}, {3,1,1,1,1}, {4,1,1,1,1}, {2,2,1,1,1},
     {3,2,1,1,1}, {4,2,1,1,1}, {3,3,1,1,1}, {4,3,1,1,1}, {4,4,1,1,1}, {2,2,2,1,1},
     {3,2,2,1,1}, {4,2,2,1,1}, {3,3,2,1,1}, {4,3,2,1,1}, {4,4,2,1,1}, {3,3,3,1,1},
     {4,3,3,1,1}, {4,4,3,1,1}, {2,2,2,2,1}, {3,2,2,2,1}, {4,2,2,2,1}, {3,3,2,2,1},
     {4,3,2,2,1}, {4,4,2,2,1}, {3,3,3,2,1}, {4,3,3,2,1}, {4,4,3,2,1}, {3,3,1,1,0}};

    int arr[] = {4,3,1,1,0},
        lower[] = {3,3,1,1,0},
        upper[] = {4,4,3,2,1};
    
    int res = 0;
    for(int i = 0; i < 48; i++) {
        assert(res == 0);
        res = _array_increment_by_limits(arr, 5, lower, upper);
        assert(memcmp(arr, test_data[i], sizeof(int) * 5) == 0);
    }
    assert(res == 1);
}

void test_gt_min_int_pattern() {
    int toprow[]  = {5,5,3,1,0},
        toprow2[] = {9,7,1,0};
    
    int result[]  = {5,5,3,1,0,5,3,1,0,3,1,0,1,0,0},
        result2[] = {9,7,1,0,7,1,0,1,0,0};
    
    int *pattern = gt_min_int_pattern(toprow, 5);
    assert(memcmp(pattern, result, 15 * sizeof(int)));
    free(pattern);

    pattern = gt_min_int_pattern(toprow2, 4);
    assert(memcmp(pattern, result2, 10 * sizeof(int)));
    free(pattern);
}

int main(int argc, char **argv) {
    test__array_increment_by_limits();
    test_gt_min_int_pattern();
}