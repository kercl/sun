from sympy.matrices import SparseMatrix

from sun.core cimport IrrepBase

cdef class Irrep(IrrepBase):
  """
  Stores the values in numpy resp. scipy arrays
  """

  def _diagonal_sparse(self, diag):
    return SparseMatrix(self.dim, self.dim, {(i, i): k for i, k in enumerate(diag)})

  def cartan(self, l):
    if self._cache_cartan[l] is None:
      diag = self._build_cartan_diagonal(l)
      self._cache_cartan[l] = self._diagonal_sparse(diag) / 2

    return self._cache_cartan[l]
