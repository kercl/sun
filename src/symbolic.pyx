from sympy.matrices import SparseMatrix
import sympy

from sun.core cimport IrrepBase

cdef class Irrep(IrrepBase):
  """
  Stores the values in numpy resp. scipy arrays
  """

  def _diagonal_sparse(self, diag):
    return SparseMatrix(self.dim, self.dim, {(i, i): k for i, k in enumerate(diag)})

  def _coord_sparse(self, entries):
    entries = {(r,c): sympy.sqrt(sympy.Rational(n,d))
              for r, c, n, d in entries}
    return SparseMatrix(self.dim, self.dim, entries)

  def cartan(self, l):
    if self._cache_cartan[l] is None:
      diag = self._build_cartan_diagonal(l)
      self._cache_cartan[l] = self._diagonal_sparse(diag) / 2

    return self._cache_cartan[l]

  def lowering(self, l):
    if self._cache_lowering[l] is None:
      entries = self._build_lowering_operator_entries(l)
      self._cache_lowering[l] = self._coord_sparse(entries)

    return self._cache_lowering[l]

  def raising(self, l):
    return self.lowering(l).adjoint()

class LieAlgebra:
  """
  TODO: Documentation
  """

  @property
  def fundamental(self):
    return self._fundamental

  def __init__(self, n):
    self._fundamental = Irrep([1] + [0] * (n - 1))