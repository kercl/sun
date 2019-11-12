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
    return coo_matrix((np.sqrt(n/d), (r, c)), shape=(self.dim, self.dim))

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
  """
  SU(3) representations
  """
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

  @staticmethod
  def gell_mann():
    irrep = SU3(1, 0)
    return [2 * irrep[i] for i in range(1, 8+1)]

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

    trx2[0] = (self._onb_csa[0].diagonal()**2).sum()
    for i in range(1, self.dim_csa):
      for j in range(i):
        trxj2 = trx2[j]
        trxixj = (self._onb_csa[i].diagonal()*
                  self._onb_csa[j].diagonal()).sum()
        
        self._onb_csa[i] -= (trxixj / trxj2) * self._onb_csa[j]

      trx2[i] = (self._onb_csa[i].diagonal()**2).sum()

    for i in range(self.dim_csa):
      self._onb_csa[i] *= np.sqrt(trx2[0] / trx2[i])

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
