from sun cimport IrrepBase

cdef class IrrepNumeric(IrrepBase):
  """
  Stores the values in numpy arrays
  """

  def cartan(self):
    if self._gt_basis.num_patterns == 0:
      self._construct_gt_basis()
