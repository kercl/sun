"""
sympy python bindings for su(n)
functionality
"""

import numpy

import sun.sun_core

class IrrepOperator:
    def __init__(self, operator, a, b):
        assert operator in ["tensor_product", "direct_sum"]

        self._factors = []

        if isinstance(a, TensorProduct):
            self._factors += a._factors
        else:
            self._factors.append(a)
        
        if isinstance(b, TensorProduct):
            self._factors += b._factors
        else:
            self._factors.append(b)

    def __str__(self):
        return u"\u2297".join(map(str, self._factors))

    def __repr__(self):
        return str(self)

    def tensor(self, repr):
        return TensorProduct(self, repr)

class Irrep:
    """
    Irreducible representation
    """

    def __init__(self, dynkin=None):
        self._dynkin = dynkin

    def __str__(self):
        """
        TODO: documentation
        """
        if self._dynkin:
            return "D(" + ",".join(map(str, self._dynkin)) + ")"

    def __repr__(self):
        """
        TODO: documentation
        """
        return str(self)

    def tensor(self, repr):
        return TensorProduct(self, repr)

    def dimension(self):
        """
        Compute the dimension of
        the irreducible representation
        """
        return sun.sun_core.irrep_dimension_from_dynkin(self._dynkin)