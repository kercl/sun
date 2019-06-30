cimport cython
from cython.view cimport array as cvarray
from libc.stdlib cimport malloc, free
from itertools import combinations

import numpy as np

cdef extern from "int_gt.h":
  cdef struct gt_tree:
    gt_int_t *array_representation
    size_t num_patterns

  size_t gt_num_of_patterns(gt_int_t *toprow, size_t length)
  void gt_free_tree(gt_tree *tree, int free_array)
  void gt_generate_all(gt_int_t **pattern, size_t *num_entries, gt_int_t *toprow, size_t length) nogil
  void gt_list_to_tree(gt_tree *tree, gt_int_t *patterns, size_t num_patterns, size_t length) nogil

cdef extern from "irrep.h":
  ctypedef int mat_int_t

  void csa_generator_diag_from_gt(gt_tree *patterns, size_t l, mat_int_t *diagonal)
  size_t lowering_operator_from_gt(gt_tree *patterns,
                                   size_t l,
                                   mat_int_t **numerators,
                                   mat_int_t **denominators,
                                   size_t **row,
                                   size_t **col)


def _l_to_pq(l):
    w = np.floor((np.sqrt(8 * l + 1) - 1) / 2)
    t = w * (w + 1) / 2
    y = l - t
    x = w - y
    return int(y + x), int(y)


def _pq_to_l(p, q):
  return int(p*(p + 1) // 2 + q)


cdef class IrrepBase:
  """
  Base class for irreducible representations
  """

  def _from_dynkin(self, dynkin):
    """
    Convert a dynkin label to the top row of the
    corresponding Gelfand-Tsetlin patterns.
    These patterns are internally used to generate
    the representations.
    """

    cdef gt_int_t *top_row

    self._length = len(dynkin) + 1
    top_row = <gt_int_t*>malloc(self._length * cython.sizeof(gt_int_t))

    if top_row is NULL:
        raise MemoryError()

    top_row[len(dynkin)] = 0
    for i in range(len(dynkin) - 1, -1, -1):
      top_row[i] = top_row[i + 1] + dynkin[i]

    self._gt_top_row = top_row

  def _from_gttoprow(self, gtoprow):
    """
    Set Gelfand-Tsetlin patterns.
    """

    cdef gt_int_t *top_row

    self._length = len(gtoprow)
    top_row = <gt_int_t*>malloc(self._length * cython.sizeof(gt_int_t))

    if top_row is NULL:
        raise MemoryError()

    for i in range(0, self._length):
      top_row[i] = gtoprow[i]

    self._gt_top_row = top_row

  cdef _construct_gt_basis(self):
    """
    construct the basis for the vector space,
    consisting of all possible Gelfand-Tsetlin
    patterns for the given top row.
    """

    cdef gt_int_t *patterns
    cdef size_t num_patterns

    with nogil:
      gt_generate_all(&patterns,
                      &num_patterns,
                      self._gt_top_row,
                      self._length)

      gt_list_to_tree(&self._gt_basis,
                      patterns,
                      num_patterns,
                      self._length)

  @property
  def dim_csa(self):
    """
    Returns the dimension of the Cartan subalgebra
    """

    return self._length - 1

  @property
  def N(self):
    return self._length
  
  @property
  def dim_lie_algebra(self):
    """
    Returns the dimension of the Lie algebra
    """

    return self._length ** 2 - 1

  @property
  def num_root_generators(self):
    """
    Returns the number of root generators
    """

    return self._length - 1

  @property
  def dim(self):
    """
    Return the dimension of the representation
    """
    if self._gt_basis.num_patterns > 0:
      return self._gt_basis.num_patterns

    cdef int d
    d = gt_num_of_patterns(self._gt_top_row, self._length)

    return int(d)

  @property
  def young(self):
    """
    return young diagram of the representation
    """

    raise NotImplementedError("TODO")

  def __init__(self, **kwargs):
    assert len(kwargs) == 1, \
      "IrrepBase takes only one of the following" \
      "keyword arguments: dynkin, gelfand_tsetlin, young"

    if "dynkin" in kwargs:
      self._from_dynkin(kwargs["dynkin"])

    self._gt_basis.num_patterns = 0

    self._cache_cartan = [None] * self.dim_csa
    self._cache_lowering = [None] * int(self.dim_csa * (self.dim_csa + 1) // 2)

    self._i = None

  def _build_cartan_diagonal(self, l):
    if self._gt_basis.num_patterns == 0:
      self._construct_gt_basis()

    cdef mat_int_t *diagonal = <mat_int_t*>malloc(self.dim * cython.sizeof(mat_int_t))
    csa_generator_diag_from_gt(&self._gt_basis, l + 1, diagonal)

    narr = np.ndarray(shape=self.dim, dtype=np.dtype("i"))
    cdef int [:] narr_view = narr

    for i in range(self.dim):
      narr_view[i] = diagonal[i]

    free(diagonal)
    return narr

  def _build_lowering_operator_entries(self, l):
    if self._gt_basis.num_patterns == 0:
      self._construct_gt_basis()

    cdef mat_int_t *numerators
    cdef mat_int_t *denominators
    cdef size_t *rows
    cdef size_t *cols
    cdef size_t num_entries

    num_entries = lowering_operator_from_gt(&self._gt_basis, l + 1,
                                            &numerators, &denominators,
                                            &rows, &cols)

    narr = np.ndarray(shape=(num_entries, 4), dtype=np.int64)
    cdef long [:, :] narr_view = narr
    for i in range(num_entries):
      narr_view[i, 0] = rows[i]
      narr_view[i, 1] = cols[i]
      narr_view[i, 2] = numerators[i]
      narr_view[i, 3] = denominators[i]

    free(numerators)
    free(denominators)
    free(rows)
    free(cols)

    return narr

  def _cartan(self, l):
    raise NotImplementedError("_cartan not implemented in IrrepBase")

  def cartan(self, l):
    if self._cache_cartan[l] is None:
      self._cache_cartan[l] = self._cartan(l)

    return self._cache_cartan[l]

  def _lowering_root(self, p):
    raise NotImplementedError("_lowering_root not implemented in IrrepBase")

  def _lowering(self, p, q):
    if p < q:
      raise KeyError(f"Unallowed configuration: p < q.")

    l = _pq_to_l(p, q)

    if p == q:
      if self._cache_lowering[l] is None:
        self._cache_lowering[l] = self._lowering_root(p)
      return self._cache_lowering[l]

    if self._cache_lowering[l] is not None:
      return self._cache_lowering[l]

    A = self._lowering(q, q)
    for k in range(q + 1, p + 1):
      B = self._lowering(k, k)
      A = B*A - A*B

    self._cache_lowering[l] = A
    return A

  def _raising(self, p, q):
    raise NotImplementedError("_raising not implemented in IrrepBase")

  def y(self, l):
    """
    Generate a basis for the the Lie algebra of
    the form: 

    Let L_i, R_i be the ith lowering and raising
    operators.

    TODO: documentation
    """

    n = self.num_root_generators
    if l >= n * (n + 1):
      return self.cartan(l - n * (n + 1))

    if l % 2 == 0:
      return self._raising(*_l_to_pq(l // 2))
    return self._lowering(*_l_to_pq((l - 1) // 2))

  def x(self, k):
    """
    Given the basis Y, compute:

    X[0] = (Y[0] + ajoint(Y[0])) / 2
    X[1] = (Y[1] - ajoint(Y[1])) / 2
    ...
    """

    n = self.num_root_generators

    if k < n * (n + 1):
      odd = k % 2

      Y1 = self.y(k - odd)
      Y2 = self.y(k - odd + 1)
      if odd == 0:
        return (Y1 + Y2) / 2
      else:
        return (Y1 - Y2) / (2 * self._i)
    else:
      return self.y(k)

  def drop_basis(self):
    if self._gt_basis.num_patterns > 0:
      gt_free_tree(&self._gt_basis, 1)
      self._gt_basis.num_patterns = 0

  def __dealloc__(self):
    self.drop_basis()
    free(self._gt_top_row)
