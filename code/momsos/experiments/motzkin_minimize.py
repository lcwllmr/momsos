import cvxpy as cp
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from momsos.poly import PolynomialRing, Polynomial
from momsos.helpers import build_plot_in_context

def motzkin_lower_bound(hierarchy_level):
    # Maximize lam such that:
    # p - lam = s0 + s1*(outer^2 - x^2 - y^2) + s2 * g1,
    # where g1 defines K

    radius = 0.2
    center = [1, -1]

    p = Polynomial.motzkin()
    polyring = PolynomialRing(p.n)
    basis0 = polyring.monomials_upto(hierarchy_level)
    basis1 = polyring.monomials_upto(hierarchy_level - 1)

    # lhs = motzkin - lam
    lam = cp.Variable()
    lhs = p
    lhs.coefficients[(0, 0)] = lhs[(0,0)] - lam

    # rhs = element of quadratic module
    Q0 = cp.Variable((len(basis0), len(basis0)), PSD=True)
    p0 = Polynomial.from_gram_matrix(p.n, hierarchy_level, Q0)
    assert p0.degree == 2 * hierarchy_level

    Q1 = cp.Variable((len(basis1), len(basis1)), PSD=True)
    p1 = Polynomial.from_gram_matrix(p.n, hierarchy_level - 1, Q1)
    assert p1.degree == 2 * (hierarchy_level - 1)
    p1 = p1 * Polynomial.ball(p.n, radius, center=center)
    assert p1.degree == 2 * hierarchy_level

    rhs = p0 + p1

    # enforce coefficient equality up to degree double hierarchy level
    constraints = []
    for m in rhs.coefficients.keys():
        constraints.append(lhs[m] == rhs[m])

    # add objective and solve
    sdp = cp.Problem(cp.Maximize(lam), constraints)
    sdp.solve()
    return sdp.status, lam.value

def main():
    ds = list(range(3, 7 + 1))
    lbs = []

    for d in ds:
        status, g = motzkin_lower_bound(d)
        lbs.append(g)
        print(f"t={d}: status={status}, lower_bound={g}")

    def plot(figure: plt.Figure):
        ax = figure.add_subplot(111)
        ax.plot(ds, lbs, marker="o")
        ax.axhline(0.0, linestyle="--", color='r', label="true minimum")
        ax.set_xlabel("hierarchy level d")
        ax.set_ylabel("lower bound $p^{\\text{sos}}_d$")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(True, alpha=0.3)
        ax.legend(loc="lower right")
    build_plot_in_context("motzkin-minimize", (7.5, 4.5), plot)

if __name__ == "__main__":
    main()

