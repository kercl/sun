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

  def __init__(self, j):
    super().__init__(dykin=1)

  def __getitem__(self, i):
    if i == 3:
      return self.cartan(0)
    return []

  @staticmethod
  def pauli():
    irrep = SU2(1)
