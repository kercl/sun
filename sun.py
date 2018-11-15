import itertools
from multiprocessing import Pool

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
  
  def unravel(self):
    flattend = self._m_raligned[::-1].flatten()
    return flattend[~np.isnan(flattend)]

class GTBasis:
  def __init__(self, max_value, min_value):
    pass
  
  @staticmethod
  def _increment_pattern_top_fixed(p):
    for k in range(len(p) - 1): # iterate through rows
      for l in range(k + 1):
        print(k,l)

  
  @classmethod
  def from_fixed_top_row(cls, row):
    pass