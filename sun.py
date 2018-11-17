import itertools
from multiprocessing import Pool

import scipy.sparse as sparse
import numpy as np
import pandas as pd

class GTPattern:
  def _init_from_array(self, pattern):
    if np.any(np.diff(np.array([len(m) for m in pattern])) != -1): raise ValueError("Pattern has incorrect shape")
    if len(pattern) != len(pattern[0]): raise ValueError("Pattern has incorrect shape")

    self._len = len(pattern)
    self._m_laligned = np.full((self._len,self._len), np.nan)
    self._m_raligned = np.full((self._len,self._len), np.nan)
    
    for i,row in enumerate(reversed(pattern)):
      self._m_laligned[i,:i+1] = row
      self._m_raligned[i,-(i+1):] = row

  def _duplicate_pattern(self, pattern):
    self._m_laligned = np.copy(pattern._m_laligned)
    self._m_raligned = np.copy(pattern._m_raligned)
    self._len = pattern._len

  def __init__(self, pattern, copy=False):
    if copy:
      self._duplicate_pattern(pattern)
    else:
      self._init_from_array(pattern)
  
  @classmethod
  def unflatten(self, array):
    N = int( (-1 + np.sqrt(1 + 8*len(array))) // 2 )

    if int(N*(N+1)//2) != len(array): raise ValueError("Invalid array. Cannot reshape into pattern.")

    indices = np.insert(np.cumsum(np.arange(N,0,-1)), 0, 0)
    return GTPattern([array[indices[i-1]:indices[i]] for i in range(1,N+1)])

  def __eq__(self, other):
    return np.equal(self._m_laligned[np.tril_indices(self._len)], other._m_laligned[np.tril_indices(other._len)]).all()

  @staticmethod
  def _align_to_left(m):
    res_m = np.full((len(m),len(m)), np.nan)
    for i in range(len(m)):
      res_m[i,:i+1] = m[i,-(i+1):]
    return res_m

  @staticmethod
  def _align_to_right(m):
    res_m = np.full((len(m),len(m)), np.nan)
    for i in range(len(m)):
      res_m[i,-(i+1):] = m[i,:i+1]
    return res_m
  
  @staticmethod
  def _tilt_left_aligned(m):
    ret = np.zeros((len(m), len(m)))
    for i,r in enumerate(m):
      ret += np.diagflat(r[:i+1], len(m) - i - 1)
    return ret

  def __getitem__(self, key):
    crop = (slice(None),slice(None))
    if not isinstance(key, tuple):
      l,k  = key,slice(None)
    else:
      l,k = key
    m = self._m_laligned
    if isinstance(k, slice) and not isinstance(l, slice):
      crop = slice(l + 1)
    elif not isinstance(k, slice) and isinstance(l, slice):
      crop = slice(k,None)
    if not isinstance(k, slice) and not isinstance(l, slice):
      return m[key]
    return m[key][crop]
  
  def __setitem__(self, key, value):
    self._m_laligned[key] = value
    
    if not np.all(np.isnan(self._m_laligned[np.triu_indices(self._len,1)])):
      self._m_laligned = GTPattern._align_to_left(self._m_raligned)
      raise IndexError("Index out of bounds.")

    return value
  
  def astype(self, t):
    ret = GTPattern(self, copy=True)
    ret._m_laligned = ret._m_laligned.astype(t)
    ret._m_raligned = ret._m_raligned.astype(t)
    return ret
  
  def _repr_html_(self):
    N = self._len
    rep = r"$\begin{pmatrix}"
    for i in range(N):
      rep += "&"*i + "&&".join(self[N - 1 - i].astype(str)) + "&"*i + r"\\"
    rep += r"\end{pmatrix}$"
    return rep

  def __str__(self):
    N = self._len
    max_digits = np.max(np.array([len(x) for x in self._m_laligned[np.tril_indices(N)].astype(str)]))
    
    rep = "gtpattern(["
    space = ""
    for i in range(N):
      rep += space + "[" + " "*max_digits*i + (" "*max_digits).join(self[N - 1 - i].astype(str)) + " "*max_digits*i + "],\n"
      space = "           "
    return rep[:-2] + "])"

  def __len__(self):
    return len(self._m_laligned)

  def sum(self, axis=None):
    r = np.nansum(self._m_laligned, axis=axis)
    if axis == 1:
      return r[::-1]
    return r
  
  def flatten(self):
    flattend = self._m_raligned[::-1].flatten()
    return flattend[~np.isnan(flattend)]

class GTBasis:
  def __init__(self, min_pattern, max_pattern):
    self._min_pattern = min_pattern
    self._max_pattern = max_pattern
  
  @staticmethod
  def _increment_array_with_limit(array, min_array, max_array):
    low_points = np.where(max_array > array)[0]

    if len(low_points) == 0:
      array[:] = np.copy(min_array)
      return array, True
    
    array[low_points[0]] += 1
    array[:low_points[0]] = array[low_points[0]]

    return array, False

  @staticmethod
  def _min_pattern_fixed_top(top_row):
    return GTPattern([top_row[i:] for i in range(len(top_row))])
  
  @staticmethod
  def _max_pattern_fixed_top(top_row):
    return GTPattern([top_row[:len(top_row)-i] for i in range(len(top_row))])

  def max(self):
    return self._max_pattern
  
  def min(self):
    return self._min_pattern
  
  def index_of(self, pattern):
    N = len(self.min())
    search_word = pattern.flatten()

    lo_bound = 0
    up_bound = len(self._pattern_map)

    for i in range(N,len(search_word)):
      lo_bound = np.searchsorted(self._pattern_map[lo_bound:up_bound,i], search_word[i]) + lo_bound
      up_bound = np.searchsorted(self._pattern_map[lo_bound:up_bound,i], search_word[i], side="right") + lo_bound

      if lo_bound + 1 == up_bound:
        return lo_bound
    
    raise IndexError("Pattern does not exist")

  def _increment_pattern_top_fixed(self, p):
    carry = True
    col_from_right = 1
    while carry:
      if col_from_right >= len(p):
        raise RuntimeError("Maximum GTPattern reached")
      
      row, col = slice(col_from_right-1,len(p)-1), -col_from_right

      new_col, carry = GTBasis._increment_array_with_limit(
        p._m_raligned[row, col],
        self.min()._m_raligned[col_from_right-1:len(p)-1,col],
        p._m_raligned[col_from_right:,col - 1]
      )
      p._m_raligned[row, col] = new_col
      if carry:
        p._m_raligned[:,-col_from_right] = self.min()._m_raligned[:,-col_from_right]
      p._m_laligned = GTPattern._align_to_left(p._m_raligned)

      col_from_right += 1
    
    return p
  
  def __len__(self):
    return self._number_of_patterns

  def iterpattern(self):
    for p in self._pattern_map:
      yield GTPattern.unflatten(p)

  @classmethod
  def from_fixed_top_row(cls, row):
    ret = cls(
      GTBasis._min_pattern_fixed_top(row),
      GTBasis._max_pattern_fixed_top(row)
    )

    ret._number_of_patterns = 1
    for k in itertools.combinations(range(len(row)),2):
      ret._number_of_patterns *= 1 + (row[k[0]] - row[k[1]]) / (k[1] - k[0])
      ret._number_of_patterns = int(ret._number_of_patterns)
    
    N = len(row)
    ret._pattern_map = np.full((ret._number_of_patterns, N*(N+1)//2), np.nan)

    p = GTPattern(ret.min(), copy=True)
    ret._pattern_map[0,:] = p.flatten()
    
    for i in range(1,ret._number_of_patterns):
      p = ret._increment_pattern_top_fixed(p)
      ret._pattern_map[i,:] = p.flatten()
    
    return ret


def matrix_scalar_prod(X,Y):
  return X.dot(Y).diagonal().sum()

class SUNLieAlgebraBasis:
  def __init__(self, rep, normalize=False, multiprocessing=True):
    self._v_basis = GTBasis.from_fixed_top_row(np.array(rep + [0]))
    self._N = len(rep) + 1
    
    self._carr_sp_dim = len(self._v_basis)

    self._build_center()

    if multiprocessing:
      p = Pool(8)
      self.raising_ops = p.map(self._build_raising_op, list(range(1,self._N)))
    else:
      self.raising_ops = []
      for l in range(1, self._N):
        self.raising_ops.append(self._build_raising_op(l))
    self._finalize_algebra()

    if normalize:
      self.onb()
  
  def casimir2(self):
    ret = 0
    for X in itertools.chain(self.center, self.center_oc):
      ret = X.dot(X) + ret
    return ret

  def _build_center(self):
    self.center = []

    row_sums = np.array([np.sum(m[:,:], axis=1) for m in self._v_basis.iterpattern()])
    sigma = np.append(row_sums,np.zeros((len(row_sums),1)),axis=1)
    center = (sigma[:,1:3] - (1/2) * (sigma[:,0:2] + sigma[:,2:4])).T
    
    self.center.append(sparse.dia_matrix((center[0],0),shape=(len(center[0]), len(center[0]))))
    self.center.append(sparse.dia_matrix((center[1],0),shape=(len(center[1]), len(center[1]))))

  def _finalize_algebra(self):
    Xall = [X for X in np.copy(self.raising_ops)]
    for X,Y in itertools.combinations(self.raising_ops, 2):
      Xall.append(X.dot(Y) - Y.dot(X))
    
    self.center_oc = []
    for X in Xall:
      self.center_oc.append(0.5 * (X + X.T))
      self.center_oc.append(0.5j * (X - X.T))
    
  def make_orthogonal(self):
    new_center = []
    scale = np.sqrt(matrix_scalar_prod(self.center[0],self.center[0]))

    for i,X in enumerate(self.center):
      for j in range(i):
        Y = self.center[j]
        X = X - matrix_scalar_prod(X,Y) * Y / matrix_scalar_prod(Y, Y)
        X = X * scale / np.sqrt(matrix_scalar_prod(X,X))
      new_center.append(X)
    self.center = new_center

  def onb(self, scale=np.sqrt(2)):
    # TODO: orthognalize, only normalize for now
    for X in self.center_oc:
      X *= scale / np.sqrt(X.dot(X).diagonal().sum())
    
    new_center = []
    for i,X in enumerate(self.center):
      for j in range(i):
        Y = self.center[j]
        X = X - matrix_scalar_prod(X,Y) * Y / matrix_scalar_prod(Y, Y)
      new_center.append(X)
    self.center = new_center
    
    for X in self.center:
      X *= scale / np.sqrt(X.dot(X).diagonal().sum())

  def _build_raising_op(self, l):
    entries = np.array([])

    for col,m in enumerate(self._v_basis.iterpattern()):
      for k in range(1,l+1):
        m_inc = GTPattern(m._m)
        m_inc._m[len(m_inc) - l, len(m_inc) - l + k - 1] += 1

        try:
          row = self._v_basis.index_of(m_inc)
        except IndexError:
          continue

        if row is None:
          continue
        
        N1 = np.prod(m[:,l+1] - m[k,l] + k - np.arange(1,l+2))
        N2 = np.prod(m[:,l-1] - m[k,l] + k - np.arange(1,l) - 1)
        D = (m[:,l] - m[k,l] + k - np.arange(1,l+1)) * (m[:,l] - m[k,l] + k - np.arange(1,l+1) - 1)
        D = np.prod(D[:k-1]) * np.prod(D[k:])

        N = N1 * N2
        if N != 0:
          entries = np.append(entries, np.array([row,col,-N/D]))
    entries = np.reshape(entries, (len(entries)//3, 3))

    return sparse.coo_matrix((np.sqrt(entries[:,2]), (entries[:,0], entries[:,1])), shape=(self._carr_sp_dim, self._carr_sp_dim))

  def _build_lowering_op(self, l):
    entries = np.array([])

    for col,m in enumerate(self._M):
      for k in range(1,l+1):
        m_inc = GTPattern(m._m)
        #M_inc._m[k+l-1,l] += 1
        m_inc._m[len(m_inc) - l, len(m_inc) - l + k - 1] -= 1

        row = self._M_idx_lookup.get(hash(m_inc))

        print(row, "\n", m_inc._m, hash(m_inc), "\n")
        
        N1 = np.prod(m[:,l+1] - m[k,l] + k - np.arange(1,l+2) + 1)
        N2 = np.prod(m[:,l-1] - m[k,l] + k - np.arange(1,l))
        D = (m[:,l] - m[k,l] + k - np.arange(1,l+1) + 1) * (m[:,l] - m[k,l] + k - np.arange(1,l+1))
        D = np.prod(D[:k]) * np.prod(D[k+1:])

        print(N1,N2,D)

        if row is None:
          continue
        
        N = N1 * N2
        if N != 0:
          entries = np.append(entries, np.array([row,col,N/D]))
    entries = np.reshape(entries, (len(entries)//3, 3))

    return sparse.coo_matrix((entries[:,2], (entries[:,0], entries[:,1])), shape=(self._carr_sp_dim, self._carr_sp_dim))
    
  def carrier_space_dim(self):
    return self._carr_sp_dim
  
  def N(self):
    return self._N
  
  def __len__(self):
    return self._N ** 2 - 1

