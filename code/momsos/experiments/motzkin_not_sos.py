import cvxpy as cp
from momsos.poly import Polynomial, PolynomialRing

def main():
    motzkin = Polynomial.motzkin()
    v_nd = PolynomialRing(motzkin.n).monomials_upto(motzkin.degree // 2)
    Q = cp.Variable((len(v_nd), len(v_nd)), PSD=True)
    pQ = Polynomial.from_gram_matrix(motzkin.n, motzkin.degree // 2, Q) # p = 

    linear_equality_constraints = []
    for m in PolynomialRing(motzkin.n).monomials_upto(motzkin.degree):
        linear_equality_constraints.append(motzkin[m] == pQ[m])

    # we don't need an optimization target because we are just interested in whether
    # there is a solution -> a feasibility problem
    sdp = cp.Problem(cp.Minimize(0), linear_equality_constraints)
    sdp.solve(solver=cp.CLARABEL, verbose=True)

if __name__ == "__main__":
    main()

