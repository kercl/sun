#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "sort.h"

#define ENABLE_TESTS
#define ABS(x) ((x) < 0 ? -x: x)

void _insert_pattern(struct gt_node *node, int *pattern, int index, int length) {
    /* patterns always start with the same length
    entries. We don't have to store them to the
    tree. */

    int *toprow = pattern; // toprow determines how many children
    int d_limit = length - 2;

    for(int d = length, i = 0; d < (length * (length + 1)) >> 1; d++) {
        if(node->children == NULL) {
            node->children = malloc(sizeof(struct gt_node) * (toprow[i] + 1));
            for(int j = 0; j <= toprow[i]; ++j) {
                node->children[j].children = NULL;
                node->children[j].index = -1;
            }
        }
        assert(pattern[d] <= pattern[i]);
        if(i == d_limit) {
            i = 0;
            d_limit--;
        }else {
            i++;
        }

        node = node->children + pattern[d] - toprow[length - 1];
    }
    node->index = index;
}

int gt_locate_in_tree(struct gt_tree *tree, int *pattern) {
    size_t length = tree->length;

    int *toprow = pattern; // toprow determines how many children
    int d_limit = length - 2;

    struct gt_node *node = &tree->root;

    for(int d = length, i = 0; d < (length * (length + 1)) >> 1; d++) {
        if(node->children == NULL || pattern[d] > pattern[i] || pattern[d] < toprow[length - 1])
            return -1;

        if(i == d_limit) {
            i = 0;
            d_limit--;
        }else {
            i++;
        }

        node = node->children + pattern[d] - toprow[length - 1];
    }
    return node->index;
}

void gt_list_to_tree(struct gt_tree *tree, int *patterns, size_t num_patterns, size_t length) {
    tree->top_row = malloc(sizeof(int) * length);
    memcpy(tree->top_row, patterns, sizeof(int) * length);
    tree->array_representation = patterns;
    tree->length = length;
    tree->num_patterns = num_patterns;
    tree->root.children = NULL;

    int m = (length * (length + 1)) >> 1;

    for(int i = 0; i < num_patterns; i++) {
        _insert_pattern(&tree->root, patterns + i * m, i, length);
    }
}

void _free_node(int *toprow, size_t d, size_t length, struct gt_node *node) {
    if(length == 0)
        return;

    if(d == length) {
        d = 0;
        length--;
    }

    for(int i; i <= toprow[d]; i++)
        _free_node(toprow, d+1, length, node->children + i);

    free(node->children);
}

void gt_free_tree(struct gt_tree *tree) {
    _free_node(tree->top_row, 0, tree->length, &tree->root);
    free(tree->array_representation);
}

void _merge(int *patterns, size_t offset, size_t n_entries, int start, int mid, int end) {
    int start_b = mid + 1, index;

    size_t cmp_n = n_entries - offset;
    int *value = malloc(sizeof(int) * cmp_n);
    
    if(memcmp(patterns + offset + n_entries * mid, 
              patterns + offset + n_entries * start_b, sizeof(int) * cmp_n) < 0)
    {
        return;
    }

    while(start <= mid && start_b <= end) {
        if(memcmp(patterns + offset + n_entries * start, 
                  patterns + offset + n_entries * start_b, sizeof(int) * cmp_n) < 0)
        { 
            start++;
        }else {
            index = start_b;

            memcpy(value, patterns + offset + n_entries * start_b, sizeof(int) * cmp_n);
            while(index != start) {
                memcpy(patterns + offset + n_entries * index, 
                       patterns + offset + n_entries * (index - 1), sizeof(int) * cmp_n);
                index--;
            }
            memcpy(patterns + offset + n_entries * start, value, sizeof(int) * cmp_n);
            
            start++; 
            mid++;
            start_b++;
        }
    }

    free(value);
}

void _merge_sort(int *patterns, size_t length, size_t n_entries, int l, int r) {
    int m;

    if(l < r) {
        m = (l + r) >> 1;

        _merge_sort(patterns, length, n_entries, l, m);
        _merge_sort(patterns, length, n_entries, m + 1, r);

        _merge(patterns, length, n_entries, l, m, r);
    }
}

void sort_pattern(int *patterns, size_t num_patterns, size_t length) {
    int n_entries = (length * (length + 1)) >> 1;

    _merge_sort(patterns, length, n_entries, 0, num_patterns - 1);
}