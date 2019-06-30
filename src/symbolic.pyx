from itertools import combinations, chain

import numpy as np
from sympy.matrices import SparseMatrix
from sympy import sqrt, Rational as frac, I

from sun.core cimport IrrepBase

cdef class Irrep(IrrepBase):
  """
  Stores the values in sympy SparseMatrix
  """

  def __init__(self, **kwargs):
    super().__init__(**kwargs)
    self._i = I

  def _diagonal_sparse(self, diag):
    return SparseMatrix(self.dim, self.dim, {(i, i): k for i, k in enumerate(diag)})

  def _coord_sparse(self, entries):
    entries = {(r,c): sqrt(frac(n,d))
              for r, c, n, d in entries}
    return SparseMatrix(self.dim, self.dim, entries)

  def _cartan(self, l):
    diag = self._build_cartan_diagonal(l)
    return self._diagonal_sparse(diag) / 2

  def _lowering_root(self, p):
    entries = self._build_lowering_operator_entries(p)
    return self._coord_sparse(entries)

  def _raising(self, p, q):
    return self._lowering(p, q).adjoint()


cdef class SU2(Irrep):
  """
  SU(2) representations
  """

  def __init__(self, j):
    super().__init__(dykin=[<int>j])

  @classmethod
  def pauli(cls):
    irrep = cls(1)


class LieAlgebra:
  """
  TODO: Documentation
  """

  @property
  def fundamental(self):
    return self._fundamental

  def __init__(self, n):
    self._fundamental = Irrep([1] + [0] * (n - 1))
