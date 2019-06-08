from scipy.sparse import dia_matrix

from sun.core cimport IrrepBase

cdef class Irrep(IrrepBase):
  """
  Stores the values in numpy resp. scipy arrays
  """

  def cartan(self, l):
    if self._cache_cartan[l] is None:
      diag = self._build_cartan_diagonal(l).astype("f") / 2.0
      self._cache_cartan[l] = dia_matrix((diag, 0), shape=(self.dim, self.dim))

    return self._cache_cartan[l]

cdef class SU2(Irrep):
  """
  SU(2) representations
  """

  def __init__(self, j):
    super().__init__(dykin=[<int>j])

  def __getitem__(self, i):
    if i == 3:
      return self.cartan(0)
    return []

  @staticmethod
  def pauli():
    irrep = SU2(1)
