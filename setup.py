from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

source_files = [
  "src/sun.pyx",
  "src/int_gt.c",
  "src/irrep.c"
]

extensions = [
  Extension("sun", 
            source_files,
            include_dirs=["."])
]

setup(
  name="sun",
  packages=["sun"],
  ext_modules=cythonize(extensions)
)
