import numpy as np
import cvxpy as cp

# Motzkin polynomial: m(x,y)=x^4 y^2 + x^2 y^4 - 3 x^2 y^2 + 1
MOTZKIN = {(4, 2): 1.0, (2, 4): 1.0, (2, 2): -3.0, (0, 0): 1.0}

def motzkin_xy(x, y):
    return x**4 * y**2 + x**2 * y**4 - 3*x**2*y**2 + 1.0

def monomial_exps_upto_deg(maxdeg):
    exps = []
    for d in range(maxdeg + 1):
        for i in range(d + 1):
            exps.append((i, d - i))  # x^i y^(d-i)
    return exps

def gram_poly_coeffs(Q, basis_exps):
    # Returns dict: exp -> cvxpy expression, for z^T Q z where z is monomial basis.
    coeffs = {}
    n = len(basis_exps)
    for i in range(n):
        ai = basis_exps[i]
        for j in range(i, n):
            aj = basis_exps[j]
            e = (ai[0] + aj[0], ai[1] + aj[1])
            term = Q[i, j] if i == j else 2 * Q[i, j]
            coeffs[e] = coeffs.get(e, 0) + term
    return coeffs

def add_scaled_shifted_poly(dst, src, scale, shift_exp):
    # dst += scale * x^shift_exp * src
    sx, sy = shift_exp
    for (ex, ey), v in src.items():
        e2 = (ex + sx, ey + sy)
        dst[e2] = dst.get(e2, 0) + scale * v
    return dst

def poly_equalities(p_expr, p_target, maxdeg):
    # Enforce coefficient-wise equality for all monomials up to total degree maxdeg.
    cons = []
    for d in range(maxdeg + 1):
        for i in range(d + 1):
            e = (i, d - i)
            lhs = p_expr.get(e, 0)
            rhs = p_target.get(e, 0.0)
            cons.append(lhs == rhs)
    return cons

def sos_decomposition_constraints(poly_target, deg_sos):
    # Find SOS s(x,y)=z^T Q z matching poly_target, with z up to degree deg_sos
    basis = monomial_exps_upto_deg(deg_sos)
    Q = cp.Variable((len(basis), len(basis)), PSD=True)
    pQ = gram_poly_coeffs(Q, basis)
    cons = poly_equalities(pQ, poly_target, maxdeg=2*deg_sos)
    return Q, cons

def solve(prob, solver="clarabel", verbose=False):
    return prob.solve(solver=cp.CLARABEL, verbose=verbose)
