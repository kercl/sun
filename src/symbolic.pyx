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

  def lowering(self, l):
    if self._cache_lowering[l] is None:
      entries = self._build_lowering_operator_entries(l)
      self._cache_lowering[l] = self._coord_sparse(entries)

    return self._cache_lowering[l]

  def raising(self, l):
    return self.lowering(l).adjoint()

  def _raising_op_from_root_gen(self, k, l):
    if k == l:
      return self.raising(k)

    A = self._raising_op_from_root_gen(k, l - 1)
    k += 1
    l -= 1
    if l < k:
      B = self.raising(k)
    else:
      B = self._raising_op_from_root_gen(k, l)

    return A*B - B*A

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

    if k < 2 * n:
      # lin com lowering/raising
      odd = k % 2
      k = int(k / 2)
      if odd == 0:
        return self.lowering(k)
      else:
        return self.raising(k)
    elif k < n * (n + 1):
      k -= 2 * n
      odd = k % 2
      comb = np.fromiter(chain.from_iterable(combinations(range(n), 2)),
                                             int, n * (n - 1))
      # round k down to highest even valued element < k
      b, a = comb[k - odd], comb[k - odd + 1]
      if odd == 0:
        return (self.lowering(a)*self.lowering(b) -
                self.lowering(b)*self.lowering(a))
      else:
        return (self.raising(a)*self.raising(b) -
                self.raising(b)*self.raising(a))
    elif k < self.dim_lie_algebra:
      k -= n * (n + 1)
      return self.cartan(k)

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