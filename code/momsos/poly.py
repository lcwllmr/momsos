from typing import Any
import math

class PolynomialRing:
    def __init__(self, n: int):
        self.n = n

    def monomials(self, d: int):
        """
        Return all monomials of degree d as exponent tuples, in
        lexicographic order. The length of the list will be equal
        to comb(n + d - 1, d) = comb(n + d - 1, n - 1).

        Lex order here means standard tuple lex order:
        (a1,...,an) > (b1,...,bn) if at the first index i where
        they differ, ai > bi.

        >>> PolynomialRing(2).monomials(3)
        [(3, 0), (2, 1), (1, 2), (0, 3)]
        >>> import math
        >>> assert len(PolynomialRing(5).monomials(10)) == math.comb(5 + 10 - 1, 10)
        """
        result = []

        def rec(i, remaining, current):
            if i == self.n - 1:
                result.append(tuple(current + [remaining]))
                return
            for k in range(remaining, -1, -1):  # descending → lex order
                rec(i + 1, remaining - k, current + [k])

        rec(0, d, [])
        return result

    def monomials_upto(self, d: int):
        """
        Returns all monomials of degree up to d - first sorted after
        degree and then lexicographically. This just concatenates
        PolynomialRing(n).monomials(d) for each d. There are exactly
        comb(n + d, n) = comb(n + d, d) such monomials.

        >>> PolynomialRing(2).monomials_upto(3)
        [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (3, 0), (2, 1), (1, 2), (0, 3)]
        >>> import math
        >>> assert len(PolynomialRing(5).monomials_upto(10)) == math.comb(10 + 5, 5)
        """
        return sum((self.monomials(i) for i in range(d+1)), start=[])
    

class Polynomial:
    def __init__(self, coefficients):
        self.coefficients = coefficients
        self.n = len(list(coefficients.keys())[0])

    def from_gram_matrix(n: int, d: int, Q: Any):
        """
        Computes polynomial coefficients from a given Gram matrix.
        For the monomial basis `v=v(n,d)` it computes `p = v.T @ Q @ v`.
        The Gram matrix `Q` must have the right dimensions (so the length of
        `v`) and must allow indexing.
        """
        len_basis = math.comb(n + d, d)
        assert Q.shape == (len_basis, len_basis), f"Q should have dimension {len_basis} for n={n} and d={d}"
        basis = PolynomialRing(n).monomials_upto(d)
        coeffs = {}
        for i, mi in enumerate(basis):
            for j, mj in enumerate(basis):
                mij = tuple(a + b for a, b in zip(mi, mj))
                coeffs[mij] = coeffs.get(mij, 0) + Q[i,j]
        return Polynomial(coeffs)

    def ball(n: int, radius: float):
        coeffs = {tuple(n * [0]): radius}
        for i in range(n):
            m = n * [0]
            m[i] = 2
            coeffs[tuple(m)] = -1.0
        return Polynomial(coeffs)

    def motzkin():
        """
        >>> import numpy as np 
        >>> x = [0.0, -1.0]
        >>> y = [0.0,  1.0]
        >>> Polynomial.motzkin()(np.vstack([x, y]))
        array([1., 0.])
        """
        return Polynomial({ (4,2): 1.0, (2,4): 1.0, (2,2): -3.0, (0,0): 1.0 })
    
    @property
    def degree(self):
        return max(sum(m) for m in self.coefficients)

    def __add__(self, other):
        """
        Add two polynomials.

        >>> p = Polynomial({(1, 0): 1, (0, 1): 1})   # x + y
        >>> q = Polynomial({(1, 0): 2, (0, 0): -1})  # 2x - 1
        >>> r = p + q
        >>> r[(1, 0)], r[(0, 1)], r[(0, 0)]
        (3, 1, -1)
        """
        assert self.n == other.n, "Polynomials must have the same number of variables"
        coeffs = dict(self.coefficients)
        for m, c in other.coefficients.items():
            coeffs[m] = coeffs.get(m, 0) + c
        return Polynomial(coeffs)

    def __mul__(self, other):
        """
        Multiply two polynomials.

        >>> p = Polynomial({(1, 0): 1, (0, 1): 1})   # x + y
        >>> q = Polynomial({(1, 0): 1})              # x
        >>> r = p * q
        >>> r[(2, 0)], r[(1, 1)]
        (1, 1)
        """
        assert self.n == other.n, "Polynomials must have the same number of variables"
        coeffs = {}
        for m1, c1 in self.coefficients.items():
            for m2, c2 in other.coefficients.items():
                m = tuple(a + b for a, b in zip(m1, m2))
                coeffs[m] = coeffs.get(m, 0) + c1 * c2
        return Polynomial(coeffs)

    def __call__(self, x):
        accum = 0
        for e, c in self.coefficients.items():
            prod = 1.0
            for i in range(self.n):
                prod *= x[i]**e[i]
            accum += c * prod
        return accum

    def __getitem__(self, m):
        return self.coefficients.get(m, 0)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
