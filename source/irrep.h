#ifndef IRREP_H
#define IRREP_H

#include "int_gt.h"
#include "sort.h"

void gt_center_generator(struct gt_tree *patterns, int l, int *numerator, int *denominator);

void gt_lowering_operator(struct gt_tree *patterns, int l);

#endif
