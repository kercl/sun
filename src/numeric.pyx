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
