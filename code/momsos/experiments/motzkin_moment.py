# NOTE: It turned out that Mosek was the most reliable solver
# for this computation. It is a commercial product but free
# for academic use. You can obtain your license at
#   https://www.mosek.com/products/academic-licenses/
# Put it into the correct directory on your system as
# indicated in the email you'll receive.

import warnings
import math
import cvxpy as cp
from momsos.poly import Polynomial, PolynomialRing

class MomentRelaxation:
    def __init__(self, n: int, d: int):
        self.n = n
        self.pr = PolynomialRing(n)
        self.d = d
        self.monomials_2d = self.pr.monomials_upto(2 * d)
        self.y = cp.Variable(len(self.monomials_2d))

    def objective(self, p: Polynomial):
        dp = math.ceil(p.degree / 2)
        assert self.d >= dp, f"Current hierarchy level {self.d} cannot capture a polynomial of degree {p.degree}"
        terms = []
        for alpha, c in p.coefficients.items():
            alphai = self.monomials_2d.index(alpha)
            terms.append(c * self.y[alphai])
        return sum(terms)

    def moment_matrix(self, g: Polynomial|None = None):
        if g is None:
            g = Polynomial.constant(self.n)
        dg = math.ceil(g.degree / 2)
        assert self.d >= dg, f"Current hierarchy level {self.d} cannot capture a polynomial of degree {g.degree}"

        local_d = self.d - dg
        assert local_d >= 0
        local_monomials = self.pr.monomials_upto(local_d)
        M = cp.Variable((len(local_monomials), len(local_monomials)), PSD=True)

        cons = []
        for i, alpha in enumerate(local_monomials):
            for j, beta in enumerate(local_monomials):
                rhs_terms = []
                for gamma, c in g.coefficients.items():
                    s = tuple(a + b + c for a, b, c in zip(alpha, beta, gamma))
                    si = self.monomials_2d.index(s)
                    rhs_terms.append(c * self.y[si])
                cons.append(M[i,j] == sum(rhs_terms))

        return M, cons

def main():
    warnings.filterwarnings("ignore", category=UserWarning, module="cvxpy")

    p = Polynomial.motzkin()
    n = p.n
    g1 = Polynomial.ball(n, 0.2, center=[1, -1])

    for d in range(6, 13 + 1):
        relaxation = MomentRelaxation(n, d)
        objective = relaxation.objective(p)
        m0, cons0 = relaxation.moment_matrix()
        m1, cons1 = relaxation.moment_matrix(g1)
        sdp = cp.Problem(cp.Minimize(objective),
                         [relaxation.y[0] == 1.0] + cons0 + cons1)

        try:
            sdp.solve(solver=cp.MOSEK)
            x_bar = [float(relaxation.y.value[1]), float(relaxation.y.value[2])]
            con_val = g1(x_bar)
            print(f"hierarchy level d = {d}")
            print(f"  -> solver status: {sdp.status}")
            print(f"  -> pmom_{d} = {sdp.value}")
            print(f"  -> xbar = {x_bar}")
            print(f"  -> constraint value g1(xbar): {con_val} ({'feasible' if con_val >= 0 else 'not feasible'})")
            p_bar = p(x_bar)
            print(f"  -> |p(xbar) - pmom_{d}| = {math.fabs(p_bar - sdp.value)}")
        except:
            print(f"hierarchy level d = {d}")
            print("  -> solver status: failed")


if __name__ == "__main__":
    main()
