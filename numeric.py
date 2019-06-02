"""
sympy python bindings for su(n)
functionality
"""

import numpy

import sun.sun_core

class Irrep:
    """
    Irreducible representation
    """

    def __init__(self, dynkin=None):
        self._dynkin = dynkin

    def __repr__(self):
        """
        TODO: documentation
        """
        if self._dynkin:
            return "D(" + ",".join(map(str, self._dynkin)) + ")"

    def dimension(self):
        """
        Compute the dimension of
        the irreducible representation
        """
        return sun.sun_core.irrep_dimension_from_dynkin(self._dynkin)