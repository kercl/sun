from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

sourcefiles = [
  "src/__init__pyx",
  "src/symbolic.pyx",
  "src/numeric.pyx",
  "int_gt.cc",
  "irrep.c"
]

extensions = [Extension("sun", sourcefiles)]

setup(
    ext_modules=cythonize(extensions)
)
