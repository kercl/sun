from distutils.core import setup, Extension

setup(name = 'sun_core', version = '1.0',  \
   ext_modules = [Extension('sun_core', [
       'sun_core/int_gt.c',
       'sun_core/irrep.c',
       'sun_core/suncoremodule.c',
   ], include_dirs=['.'])])