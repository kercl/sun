#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "int_gt.h"

int _array_increment_by_limits(int *arr, size_t length, int *lower_lim, int *upper_lim);
int _gt_increment_transposed(int *pattern_tr, int *min_pattern_tr, size_t length);

void print_pattern_raligned(int *pattern, size_t length) {
    for(int i = 0; i < length; ++i) {
        for(int k = 0; k < i; k++)
            printf("   ");
        for(int j = i; j < length; ++j) {
            printf("%2d ", pattern[GT_RL_IDX(i, j, length)]);
        }
        printf("\n");
    }
}

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

void test_gt_allocate_min_int_pattern() {
    int toprow[]  = {5,5,3,1,0},
        toprow2[] = {9,7,1,0};
    
    int result[]  = {5,5,3,1,0,5,3,1,0,3,1,0,1,0,0},
        result2[] = {9,7,1,0,7,1,0,1,0,0};
    
    int *pattern = gt_allocate_min_int_pattern(toprow, 5);
    assert(memcmp(pattern, result, 15 * sizeof(int)) == 0);
    free(pattern);

    pattern = gt_allocate_min_int_pattern(toprow2, 4);
    assert(memcmp(pattern, result2, 10 * sizeof(int)) == 0);
    free(pattern);
}

void test_gt_transpose() {
    int tr_1[] = {6,5,4,3,2,1,0},
        tr_2[] = {8,7,6,5,4};

    int res_1[] = {
         0, 0, 0, 0, 0, 0, 0, 
            1, 1, 1, 1, 1, 1,
               2, 2, 2, 2, 2, 
                  3, 3, 3, 3, 
                     4, 4, 4, 
                        5, 5, 
                           6};
    int res_2[] = {
         4, 4, 4, 4, 4, 
            5, 5, 5, 5, 
               6, 6, 6, 
                  7, 7, 
                     8};
        
    int *pattern = gt_allocate_min_int_pattern(tr_1, 7);
    gt_transpose(pattern, 7);
    assert(memcmp(pattern, res_1, 28*sizeof(int)) == 0);
    free(pattern);

    pattern = gt_allocate_min_int_pattern(tr_2, 5);
    gt_transpose(pattern, 5);
    assert(memcmp(pattern, res_2, 15*sizeof(int)) == 0);
    free(pattern);
}

void test_gt_generate_all_transposed() {
    int toprow[] = {300,0,0}, length = 3;
    
    int *patterns,
         n_entries;

    gt_generate_all_transposed(&patterns, &n_entries, toprow, length);
    printf("%d\n", n_entries);
    // int m = length * (length + 1) >> 1;
    // for(int i = 0; i < n_entries; ++i)
    //     print_pattern_raligned(patterns + m * i, length);
}

int main(int argc, char **argv) {
    test__array_increment_by_limits();
    test_gt_allocate_min_int_pattern();
    test_gt_transpose();

    int toprow[] = {2,1,0};
    int length = 3;
    /*int *pattern = gt_allocate_min_int_pattern(toprow, length);
    int *min_pattern = gt_allocate_min_int_pattern(toprow, length);
    
    gt_transpose(min_pattern, length);
    gt_transpose(pattern, length);

    print_pattern_raligned(pattern, length);
    printf("---\n");
    size_t dim = 1;
    for(;;) {
        int inc = _gt_increment_transposed(pattern, min_pattern, length);
        if(inc == 0)
            break;
        dim++;
    }
    print_pattern_raligned(pattern, length);

    printf("dim(V)=%ld\n", dim);*/
    printf("dim(V) = %ld\n", gt_num_of_patterns(toprow, length));
    test_gt_generate_all_transposed();
}
