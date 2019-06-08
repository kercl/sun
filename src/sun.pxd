cdef extern from "int_gt.h":
  ctypedef int gt_int_t

  cdef struct gt_tree:
    size_t num_patterns

cdef class IrrepBase:
    cdef gt_int_t *_gt_top_row
    cdef size_t _length
    cdef gt_tree _gt_basis
    cdef list _cache_cartan

    cdef _construct_gt_basis(self)
