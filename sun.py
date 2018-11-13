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
  
  def _inc_at_coords_(self, coords):
    self._m_laligned[coords] += 1
    self._m_raligned[coords[0], self._len - 1 - coords[0] + coords[1]] += 1

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

  def _get_inc_candidates(self):
    lu_inc_possible = np.vstack([np.diff(self._m_laligned, axis=0), np.full(self._len,np.nan)])

    d_inc_possible = np.vstack([np.diff(self._m_raligned[::-1], axis=0), np.full(self._len,np.nan)])
    np.fill_diagonal(d_inc_possible,1)
    d_inc_possible = GTPattern._align_to_left(d_inc_possible[::-1])

    lu_inc_possible[-1,:] = 0 # top row not icrementable

    lu_inc_possible = np.nan_to_num(lu_inc_possible)
    d_inc_possible = np.nan_to_num(d_inc_possible)

    coords = np.array(np.where(np.greater(lu_inc_possible, 0) & np.greater(d_inc_possible, 0))).T
    return [(x,y) for x,y in coords]

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

  def sum(self, axis=None):
    r = np.nansum(self._m_laligned, axis=axis)
    if axis == 1:
      return r[::-1]
    return r
  
  def increment(self, keep_top_row=False):
    print(np.diff(self._m_laligned))

class GTPatternIterator:
  MULTIPROCESSING_LIMIT = 5
  POOLSIZE = 8

  def __init__(self, start, multiprocessing=False):
    self._base_patterns = [start]
    self._current_pos = 0
    self._multiprocessing = multiprocessing
  
  def __iter__(self):
    return self

  @staticmethod
  def _digrams_for_next_weight(old_pattern):
    inc_coords = old_pattern._get_inc_candidates()
    pattern_batch = np.empty(len(inc_coords), dtype=object)
    for i,c in enumerate(inc_coords):
      pattern_batch[i] = GTPattern(old_pattern, copy=True)
      pattern_batch[i]._inc_at_coords_(c)
    return pattern_batch

  def __next__(self):
    if self._current_pos < len(self._base_patterns):
      m = self._base_patterns[self._current_pos]
      self._current_pos += 1
      return m
    
    if self._current_pos < GTPatternIterator.MULTIPROCESSING_LIMIT or not self._multiprocessing:
      new_patterns = []
      
      pattern_batch = None
      for old_pattern in self._base_patterns:
        inc_coords = old_pattern._get_inc_candidates()
        pattern_batch = np.empty(len(inc_coords), dtype=object)
        for i,c in enumerate(inc_coords):
          pattern_batch[i] = GTPattern(old_pattern, copy=True)
          pattern_batch[i]._inc_at_coords_(c)
        new_patterns.append(pattern_batch)
    else:
      with Pool(GTPatternIterator.POOLSIZE) as p:
        new_patterns = p.map(GTPatternIterator._digrams_for_next_weight, self._base_patterns)
    
    self._base_patterns = []
    for p in itertools.chain(*new_patterns):
      exists = False
      for q in self._base_patterns:
        if p == q:
          exists = True
          break
      if not exists:
        self._base_patterns.append(p)

    if len(self._base_patterns) == 0:
      raise StopIteration()
    self._current_pos = 1
    return self._base_patterns[0]