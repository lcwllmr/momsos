import math
import cvxpy as cp
from momsos.poly import Polynomial, PolynomialRing

def attempt_sos_decomposition(target):
    d1 = math.ceil(target.degree / 2)
    d2 = 2 * d1
    v_nd = PolynomialRing(target.n).monomials_upto(d1)
    Q = cp.Variable((len(v_nd), len(v_nd)), PSD=True)
    pQ = Polynomial.from_gram_matrix(target.n, d1, Q)

    linear_equality_constraints = []
    for m in PolynomialRing(target.n).monomials_upto(d2):
        linear_equality_constraints.append(target[m] == pQ[m])

    # NOTE: we don't need an optimization target because just need
    # to know whether there is a solution (feasibility problem)
    sdp = cp.Problem(cp.Minimize(0), linear_equality_constraints)
    sdp.solve(solver=cp.CLARABEL)

    print("  ...  Computation finished.")
    print("  Solution status:", sdp.status)
    print(f"  Solve time: {sdp.solver_stats.solve_time}s")
    print(f"  Iterations: {sdp.solver_stats.num_iters}")
    

def main():
    print('\nStarting test: Is M(x,y) a sum of squares?', end='')
    motzkin = Polynomial.motzkin()
    attempt_sos_decomposition(motzkin)

    print('\nStarting test: Is (x^2 + y^2) * M(x,y) a sum of squares?', end='')
    weighted_motzkin = Polynomial({(2,0): 1.0, (0,2): 1.0}) * Polynomial.motzkin()
    attempt_sos_decomposition(weighted_motzkin)

if __name__ == "__main__":
    main()

