from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

sourcefiles = [
  "src/sun.pyx",
  "src/symbolic.pyx",
  "src/numeric.pyx",
  "src/int_gt.c",
  "src/irrep.c"
]

extensions = [Extension("sun", sourcefiles)]

setup(
    ext_modules=cythonize(extensions)
)
