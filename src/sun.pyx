cimport cython
from libc.stdlib cimport malloc, free

cdef extern from "int_gt.h":
  ctypedef short int gt_int_t
  size_t gt_num_of_patterns(gt_int_t *toprow, size_t length)

cdef class IrrepBase:
  """
  Base class for irreducible representations
  """

  cdef gt_int_t *_gt_top_row
  cdef size_t _length

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

  @property
  def dim(self):
    """
    Return the dimension of the representation
    """

    return int(gt_num_of_patterns(self._gt_top_row, self._length))

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

  def __dealloc__(self):
    free(self._gt_top_row)
