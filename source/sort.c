#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "sort.h"

void _merge(int *patterns, size_t n_entries, int start, int mid, int end) {
    int start_b = mid + 1, index;
    int *value = malloc(sizeof(int) * n_entries);
    
    if(memcmp(patterns + n_entries * mid, 
              patterns + n_entries * start_b, sizeof(int) * n_entries) < 0)
    {
        return;
    }

    while(start <= mid && start_b <= end) {
        if(memcmp(patterns + n_entries * start, 
                  patterns + n_entries * start_b, sizeof(int) * n_entries) < 0)
        { 
            start++;
        }else {
            index = start_b;

            memcpy(value, patterns + n_entries * start_b, sizeof(int) * n_entries);
            while(index != start) {
                memcpy(patterns + n_entries * index, 
                       patterns + n_entries * (index - 1), sizeof(int) * n_entries);
                index--;
            }
            memcpy(patterns + n_entries * start, value, sizeof(int) * n_entries);
            
            start++; 
            mid++;
            start_b++;
        }
    }

    free(value);
}

void _merge_sort(int *patterns, size_t n_entries, int l, int r) {
    int m;

    if(l < r) {
        m = (l + r) >> 1;

        _merge_sort(patterns, n_entries, l, m);
        _merge_sort(patterns, n_entries, m + 1, r);

        _merge(patterns, n_entries, l, m, r);
    }
}

void sort_pattern(int *patterns, size_t num_patterns, size_t length) {
    int n_entries = (length * (length + 1)) >> 1;

    _merge_sort(patterns, n_entries, 0, num_patterns - 1);
}