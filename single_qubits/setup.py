from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(["hello_world.pyx","numpy_test.pyx","qutip_diag_test.pyx"],annotate=True)
)
