cimport cython
from libc.stdlib cimport malloc, free

cdef extern from "int_gt.h":
  cdef struct gt_tree:
    size_t num_patterns

  size_t gt_num_of_patterns(gt_int_t *toprow, size_t length) nogil
  void gt_generate_all(gt_int_t **pattern, size_t *num_entries, gt_int_t *toprow, size_t length) nogil


cdef extern from "irrep.h":
  ctypedef int mat_int_t

  void csa_generator_diag_from_gt(gt_tree *patterns, size_t l, mat_int_t *diagonal)


cdef class CartanSubalgebra:
  """
  Generates the Cartan subalgebra for a given Irrep
  """

  def __getitem__(self, i):
    if self._cache[i] is not None:
      return self._cache[i]

  def __init__(self, irrep):
    self._irrep = irrep
    self._cache = [None] * irrep.dim_csa


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

    for i in range(6):
      print(patterns[i])

  @property
  def dim_csa(self):
    """
    Returns the dimension of the Cartan subalgebra
    """

    return self._length - 1

  @property
  def dim(self):
    """
    Return the dimension of the representation
    """

    cdef int d
    with nogil:
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

  def __dealloc__(self):
    free(self._gt_top_row)

cdef class IrrepNumeric(IrrepBase):
  """
  Stores the values in numpy arrays
  """

  def cartan(self):
    if self._gt_basis.num_patterns == 0:
      self._construct_gt_basis()
