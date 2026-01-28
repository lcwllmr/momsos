import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from helper import MOTZKIN, monomial_exps_upto_deg, gram_poly_coeffs, add_scaled_shifted_poly, poly_equalities, solve

def lower_bound_on_ball(level_d, R2):
    # Maximize gamma such that:
    # m - gamma = s0 + s1*(R2 - x^2 - y^2),
    # with s0 SOS of deg<=2d and s1 SOS of deg<=2(d-1).
    basis0 = monomial_exps_upto_deg(level_d)
    basis1 = monomial_exps_upto_deg(level_d - 1)

    Q0 = cp.Variable((len(basis0), len(basis0)), PSD=True)
    Q1 = cp.Variable((len(basis1), len(basis1)), PSD=True)
    gamma = cp.Variable()

    p0 = gram_poly_coeffs(Q0, basis0)        # deg <= 2d
    p1 = gram_poly_coeffs(Q1, basis1)        # deg <= 2(d-1)

    # Build p1*(R2 - x^2 - y^2)
    p1g = {}
    add_scaled_shifted_poly(p1g, p1, R2, (0, 0))
    add_scaled_shifted_poly(p1g, p1, -1.0, (2, 0))
    add_scaled_shifted_poly(p1g, p1, -1.0, (0, 2))

    # Target: MOTZKIN - gamma
    target = dict(MOTZKIN)
    target[(0, 0)] = target.get((0, 0), 0.0) - gamma

    # Enforce coefficient equality up to degree 2d
    lhs = dict(p0)
    for e, v in p1g.items():
        lhs[e] = lhs.get(e, 0) + v
    cons = poly_equalities(lhs, target, maxdeg=2*level_d)

    prob = cp.Problem(cp.Maximize(gamma), cons)
    solve(prob, solver="scs", verbose=False)
    return prob.status, gamma.value

def main():
    R2 = 2  # radius^2 = 2 includes the true minimizers (±1,±1); true min is 0
    ds = list(range(2, 8))  # small but illustrative
    lbs = []

    for d in ds:
        status, g = lower_bound_on_ball(d, R2)
        if status not in ("optimal", "optimal_inaccurate"):
            g = np.nan
        lbs.append(g)
        print(f"d={d}: status={status}, lower_bound={g}")

    fig = plt.figure(figsize=(7.5, 4.5))
    ax = fig.add_subplot(111)
    ax.plot(ds, lbs, marker="o")
    ax.axhline(0.0, linestyle="--")
    ax.set_title("SOS hierarchy lower bounds for Motzkin on x^2+y^2 ≤ 2")
    ax.set_xlabel("hierarchy level d")
    ax.set_ylabel("lower bound γ_d")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("motzkin_ball_hierarchy.png", dpi=200)
    print("Wrote motzkin_ball_hierarchy.png")

if __name__ == "__main__":
    main()

