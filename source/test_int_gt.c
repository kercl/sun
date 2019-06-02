#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "int_gt.h"
#include "sort.h"
#include "irrep.h"

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

void print_pattern_flattened(int *pattern, size_t length) {
    for(int i = 0; i < (length*(length + 1)) >> 1; ++i) {
        printf("%2d ", pattern[i]);
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
    int toprow[] = {3,2,1,0};
    size_t length = 4;
    
    int *patterns;
    size_t n_entries;
    
    printf("generating patterns... (%ld)\n", gt_num_of_patterns(toprow, length));
    gt_generate_all(&patterns, &n_entries, toprow, length);
    printf("generating search tree...\n");
    struct gt_tree tree;
    gt_list_to_tree(&tree, patterns, n_entries, length);
    
    int m = length * (length + 1) >> 1;
    
    printf("testing search tree...\n");
    for(int i = 0; i < n_entries; ++i) {
        assert(i == gt_locate_in_tree(&tree, patterns + i * m));
    }

    printf("number of patterns %ld\n", n_entries);

    gt_free_tree(&tree);
}

void test_sort() {
    int patterns[] = {
        5,3,1,
        5,1,3,
        5,2,1,
        3,6,1,
        2,3,1,
        6,4,1,
        1,1,1
    };
    int length = 2;
    int n_patterns = 7;
    int m = (length * (length + 1)) >> 1;

    gt_sort_patterns(patterns, n_patterns, length);

    for(int i = 0; i < n_patterns; ++i) {
        print_pattern_flattened(patterns + m * i, length);
        printf("\n");
    }
}

int main(int argc, char **argv) {
    //test__array_increment_by_limits();
    //test_gt_allocate_min_int_pattern();
    //test_gt_transpose();
    //test_gt_generate_all_transposed();

    int toprow[] = {2,0,0};
    int length = sizeof(toprow) / sizeof(int);

    int *patterns;
    size_t n_patterns;

    gt_generate_all(&patterns, &n_patterns, toprow, length);
    struct gt_tree tree;
    gt_list_to_tree(&tree, patterns, n_patterns, length);

    printf("%ld\n", n_patterns);

    for(int i = 0; i < n_patterns; i++)
        print_pattern_raligned(patterns + i * (length * (length + 1) / 2), length);

    int *x3 = malloc(sizeof(int) * n_patterns * 2);
    
    gt_center_generator(&tree, 1, x3, x3 + n_patterns);

    for(int i = 0; i < n_patterns; i++) {
        printf("%d/%d ", x3[i], x3[n_patterns + i]);
    }
    printf("\n");

    int *a = malloc(sizeof(int) * n_patterns * 2);

    printf("--\n");
    gt_lowering_operator(&tree, 2);
    printf("--\n");
    gt_lowering_operator(&tree, 1);

    //printf("dim(V) = %ld\n", gt_num_of_patterns(toprow, length));
    //test_gt_generate_all_transposed();
    printf("\n");
}
