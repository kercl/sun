from scipy.sparse import dia_matrix, coo_matrix
import numpy as np

from sun.core cimport IrrepBase

cdef class Irrep(IrrepBase):
  """
  Stores the values in numpy resp. scipy arrays
  """

  def __init__(self, **kwargs):
    super().__init__(**kwargs)
    self._i = 1j

  def _cartan(self, l):
    diag = self._build_cartan_diagonal(l).astype("f") / 2.0
    return dia_matrix((diag, 0), shape=(self.dim, self.dim))

  def _lowering_root(self, p):
    r, c, n, d = self._build_lowering_operator_entries(p).T
    return coo_matrix((n/d, (r, c)), shape=(self.dim, self.dim))

  def _raising(self, p, q):
    return self._lowering(p, q).H


cdef class SU2(Irrep):
  """
  SU(2) representations
  """

  def __init__(self, n):
    super().__init__(dynkin=[n])

  def __getitem__(self, i):
    i -= 1
    return self.x(i)

  @staticmethod
  def pauli():
    irrep = SU2(1)
    return [2 * irrep.x(0),
            2 * irrep.x(1),
            2 * irrep.x(2)]

cdef class SU3(Irrep):
  cdef object _x_8

  def __init__(self, n, m):
    super().__init__(dynkin=[n, m])

    x_3 = self.x(6)
    x_8 = self.x(7)

    # optimized Tr(X_3.X_8)
    # both are diagonal matrices
    x3_norm = (x_3.diagonal() * x_3.diagonal()).sum()

    # gram-schmidt
    factor = (x_3.diagonal() * x_8.diagonal()).sum() / x3_norm
    self._x_8 = x_8 - factor * x_3

    x8_norm = (self._x_8.diagonal() * self._x_8.diagonal()).sum()

    # rescale to same norm of x_3
    self._x_8 = np.sqrt(x3_norm / x8_norm) * self._x_8

  def __getitem__(self, i):
    mapping = {
      1: 0,
      2: 1,
      3: 6,
      4: 2,
      5: 3,
      6: 4,
      7: 5
    }
    if i in mapping:
      return self.x(mapping[i])

    if i == 8:
      return self._x_8
    
    raise KeyError("Index out of bounds.")
