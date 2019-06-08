from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

source_files_core = [
  "src/core.pyx",
  "src/int_gt.c",
  "src/irrep.c"
]

source_files_numeric = [
  "src/numeric.pyx"
]

source_files_symbolic = [
  "src/symbolic.pyx"
]

extensions = [
  Extension("sun.core", 
            source_files_core,
            include_dirs=["."]),
  Extension("sun.numeric", 
            source_files_numeric,
            include_dirs=["."]),
  Extension("sun.symbolic", 
            source_files_symbolic,
            include_dirs=["."])
]

setup(
  name="sun",
  packages=["sun"],
  ext_modules=cythonize(extensions)
)
