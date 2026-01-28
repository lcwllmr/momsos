import cvxpy as cp
from helper import MOTZKIN, sos_decomposition_constraints, solve

def main():
    # If Motzkin were SOS, we'd have m(x,y)=z^T Q z with z monomials up to degree 3.
    Q, cons = sos_decomposition_constraints(MOTZKIN, deg_sos=3)
    prob = cp.Problem(cp.Minimize(0), cons)

    solve(prob, verbose=False)

    print("Motzkin SOS feasibility (degree-3 Gram) via Clarabel")
    print("status:", prob.status)
    if prob.status in ("infeasible", "infeasible_inaccurate"):
        print("Conclusion: infeasible -> Motzkin is not SOS (over R[x,y]).")
    else:
        print("Unexpected status. If 'optimal', you found an SOS Gram matrix (should not happen).")

if __name__ == "__main__":
    main()

