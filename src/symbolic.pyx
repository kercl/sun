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

cdef class LieAlgebra(Irrep):
  """
  SU(3) representations
  """
  cdef object _onb_csa

  def __init__(self, *args):
    super().__init__(dynkin=args)

    self._onb_csa = [self.cartan(i) for i in range(self.dim_csa)]
    trx2 = [None] * self.dim_csa

    # use fact that in the generated bases
    # CSA is given in diagonal form
    # tr(X_j.X_j) = sum((X_j)_kk * (X_j)_kk)

    trx2[0] = sum([self._onb_csa[0][k,k]**2 for k in range(self.dim)])
    for i in range(1, self.dim_csa):
      for j in range(i):
        trxj2 = trx2[j]
        trxixj = sum([self._onb_csa[i][k,k]*
                      self._onb_csa[j][k,k] for k in range(self.dim)])
        
        self._onb_csa[i] -= (trxixj / trxj2) * self._onb_csa[j]

      trx2[i] = sum([self._onb_csa[i][k,k]**2 for k in range(self.dim)])

    for i in range(self.dim_csa):
      self._onb_csa[i] *= sqrt(trx2[0] / trx2[i])

  @staticmethod
  def _int_sqrt(n):
    lbound = 0
    ubound = n
    
    while lbound != ubound:
      mid = (lbound + ubound) // 2
      mid2 = mid**2
      
      if mid2 == n:
        return mid, 0
      elif mid2 < n:
        lbound = mid + 1
      else:
        ubound = mid

    m = lbound**2
    if m == n:
      return lbound, 0
    return lbound - 1, n - m + 2*lbound - 1

  def __getitem__(self, i):
    if i < 1 or i > self.dim_lie_algebra:
      raise KeyError("Index out of bounds.")

    j, d = LieAlgebra._int_sqrt(i + 1)

    # i + 1 = j^2 => element is in CSA
    if d == 0:
      return self._onb_csa[j - 2]

    return self.x(i - j)
