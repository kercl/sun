#ifndef SORT_H
#define SORT_H

#include <stddef.h>

struct gt_node {
    struct gt_node *children;
    int index;
};

struct gt_tree {
    int *array_representation,
        *top_row;
    struct gt_node root;
    size_t length, num_patterns;
};

/*
 *
 */
void gt_list_to_tree(struct gt_tree *tree, int *patterns, size_t num_patterns, size_t length);

/*
 *
 */
void gt_free_tree(struct gt_tree *tree);

/*
 *
 */
int gt_locate_in_tree(struct gt_tree *tree, int *pattern);

/*
 *
 */
void sort_pattern(int *patterns, size_t num_patterns, size_t length);

#endif
