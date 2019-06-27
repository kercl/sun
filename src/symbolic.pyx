from itertools import combinations, chain

import numpy as np
from sympy.matrices import SparseMatrix
from sympy import sqrt, Rational as frac, I

from sun.core cimport IrrepBase

cdef class Irrep(IrrepBase):
  """
  Stores the values in numpy resp. scipy arrays
  """

  def _diagonal_sparse(self, diag):
    return SparseMatrix(self.dim, self.dim, {(i, i): k for i, k in enumerate(diag)})

  def _coord_sparse(self, entries):
    entries = {(r,c): sqrt(frac(n,d))
              for r, c, n, d in entries}
    return SparseMatrix(self.dim, self.dim, entries)

  def cartan(self, l):
    if self._cache_cartan[l] is None:
      diag = self._build_cartan_diagonal(l)
      self._cache_cartan[l] = self._diagonal_sparse(diag) / 2

    return self._cache_cartan[l]

  def _root_lowering(self, l):
    if self._cache_lowering[l] is None:
      entries = self._build_lowering_operator_entries(l)
      self._cache_lowering[l] = self._coord_sparse(entries)

    return self._cache_lowering[l]

  def _root_raising(self, l):
    return self.lowering(l).adjoint()

  def _lowering(self, x, y, func):
    if self._cache_lowering[l] is not None:
      return self._cache_lowering[l]

    p = int(y + (x + y) * (x + y + 1) / 2)
    
    root_gen_idx_start = (self.dim_lie_algebra -
                          self.num_root_generators)

    if p > root_gen_idx_start:
      p -= root_gen_idx_start
      self._cache_lowering[p] = self._build_lowering_operator_entries(p)
    else:
      A = self.

      A*B - B*A

  def Y(self, k):
    """
    Generate a basis for the the Lie algebra of
    the form:

    Let L_i, R_i be the ith lowering and raising
    operators, H_i the ith cartan generator and
    N the number of root generators.
    Then:

    Y_0 = L_0
    Y_1 = R_0
    ...
    Y_{2N-2} = L_{N-1}
    Y_{2N-1} = R_{N-1}
    Y_{2N} = [L_0, L_1]
    Y_{2N} = [R_0, R_1]
    ...
    Y_{N(N-1)-2} = [L_{N-2}, L_{N-1}]
    Y_{N(N-1)-1} = [R_{N-2}, R_{N-1}]
    Y_{N(N-1)} = H_0
    ...
    Y_{N^2-2} = H_{N-1}
    """
    n = self.num_root_generators
    comb = np.fromiter(chain.from_iterable(combinations(range(n), 2)),
                                           int, n * (n - 1))

    if k < 2 * n:
      odd = k % 2
      k = int(k / 2)
      if odd == 0:
        return self.raising(k)
      else:
        return self.lowering(k)
    
    k -= 2 * n
    if k < len(comb):
      odd = k % 2
      k = int(k / 2)
      if odd == 0:
        return self._lr_from_root_gen(*comb[k], func=self.raising)
      else:
        return self.lowering(k)


  def X(self, k):
    """
    Given the basis Y, compute:

    X[0] = (Y[0] + ajoint(Y[0])) / 2
    X[1] = (Y[1] - ajoint(Y[1])) / 2
    ...
    """

    n = self.num_root_generators

    if k < n * (n + 1):
      odd = k % 2

      Y = self.Y(k - odd)
      if odd == 0:
        return (Y + Y.adjoint()) / 2
      else:
        return (Y - Y.adjoint()) / (2 * I)
    else:
      return self.Y(k)

class LieAlgebra:
  """
  TODO: Documentation
  """

  @property
  def fundamental(self):
    return self._fundamental

  def __init__(self, n):
    self._fundamental = Irrep([1] + [0] * (n - 1))