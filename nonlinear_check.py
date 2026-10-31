import argparse, json, time
import numpy as np
try:
    from pypardiso import spsolve
    fmt = "csr"
except ImportError:
    from scipy.sparse.linalg import spsolve
    fmt = "csc"
from skfem import MeshHex, Basis, ElementHex2, ElementVector, BilinearForm, LinearForm, FacetBasis, solve, condense

ap = argparse.ArgumentParser()
ap.add_argument("--side", type=float, default=600.0); ap.add_argument("--t", type=float, default=13.0)
ap.add_argument("--p", type=float, default=100.0); ap.add_argument("--nxy", type=int, default=12)
ap.add_argument("--nz", type=int, default=2); ap.add_argument("--steps", type=int, default=4)
ap.add_argument("--pi", type=float, nargs=3, default=[6.6, -1.1, 138.1])
ap.add_argument("--resL", type=float, nargs=3, default=[30, 5, 5]); ap.add_argument("--resT", type=float, nargs=3, default=[30, 5, 5])
args = ap.parse_args()
t0 = time.time(); um = 1e-6
a, T, P = args.side / 2 * um, args.t * um, args.p * 1e3

def cubic(c11, c12, c44):
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if i == j and k == l: C[i, j, k, l] = c11 if i == k else c12
                    elif ((i == k and j == l) or (i == l and j == k)) and i != j: C[i, j, k, l] = c44
    return C
c, s_ = np.cos(np.pi / 4), np.sin(np.pi / 4)
R = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])
Ct = np.einsum('ia,jb,kc,ld,abcd->ijkl', R, R, R, R, cubic(165.7e9, 63.9e9, 79.6e9))
NZ = [(i, j, k, l, Ct[i, j, k, l]) for i in range(3) for j in range(3) for k in range(3) for l in range(3)
      if abs(Ct[i, j, k, l]) > 1e-6 * abs(Ct).max()]

xs = a * np.sin(np.pi / 2 * np.linspace(0, 1, args.nxy + 1))
m = MeshHex.init_tensor(xs, xs, np.linspace(0, T, args.nz + 1))
e = ElementVector(ElementHex2()); basis = Basis(m, e)
I3 = np.eye(3)

def kin(G):
    F = G + I3[:, :, None, None]
    E = 0.5 * (np.einsum('ki...,kj...->ij...', F, F) - I3[:, :, None, None])
    S = np.zeros_like(E)
    for i, j, k, l, cc in NZ: S[i, j] += cc * E[k, l]
    return F, S

def CC(X):
    out = np.zeros_like(X)
    for i, j, k, l, cc in NZ: out[i, j] += cc * X[k, l]
    return out

@BilinearForm
def tangent(du, v, w):
    F, S = w['F'], w['S']
    gdu, gv = du.grad, v.grad
    dE = 0.5 * (np.einsum('ki...,kj...->ij...', gdu, F) + np.einsum('ki...,kj...->ij...', F, gdu))
    dS = CC(dE)
    dP = np.einsum('ik...,kj...->ij...', gdu, S) + np.einsum('ik...,kj...->ij...', F, dS)
    return np.einsum('ij...,ij...->...', dP, gv)

@LinearForm
def internal(v, w):
    Pk = np.einsum('ik...,kj...->ij...', w['F'], w['S'])
    return np.einsum('ij...,ij...->...', Pk, v.grad)

@LinearForm
def load(v, w):
    return -v.value[2]

fb = FacetBasis(m, e, facets=m.facets_satisfying(lambda x: np.isclose(x[2], T)))
f1 = load.assemble(fb)
D = basis.get_dofs
fixed = np.unique(np.concatenate([D(m.facets_satisfying(lambda x: np.isclose(x[0], 0))).all('u^1'),
                                  D(m.facets_satisfying(lambda x: np.isclose(x[1], 0))).all('u^2'),
                                  D(m.facets_satisfying(lambda x: np.isclose(x[0], a) | np.isclose(x[1], a))).all()]))
def lin(A, b):
    return spsolve(A.asformat(fmt), b)

u = np.zeros(basis.N); log = []
for step in range(1, args.steps + 1):
    load_now = P * step / args.steps
    for it in range(25):
        F, S = kin(basis.interpolate(u).grad)
        r = load_now * f1 - internal.assemble(basis, F=F, S=S)
        free = np.setdiff1d(np.arange(basis.N), fixed)
        res = np.linalg.norm(r[free]) / np.linalg.norm(load_now * f1[free])
        if res < 1e-8: break
        K = tangent.assemble(basis, F=F, S=S)
        u += solve(*condense(K, r, D=fixed), solver=lin)
    log.append((step, it, float(res)))

cd = basis.get_dofs(nodes=lambda x: np.isclose(x[0], 0) & np.isclose(x[1], 0)).nodal['u^3']
w0 = float(np.mean(np.abs(u[cd])))
F, S = kin(basis.interpolate(u).grad)
sb = basis.with_element(ElementHex2())
@BilinearForm
def mass(p_, q_, w): return p_ * q_
@LinearForm
def rhs(q_, w): return w['s'] * q_
M = mass.assemble(sb)
nod = {n: lin(M, rhs.assemble(sb, s=S[i, j])) for (i, j), n in {(0, 0): 'sxx', (1, 1): 'syy'}.items()}
def avg(x0, x1, y0, y1, n=25):
    X, Y = np.meshgrid(np.linspace(x0, x1, n), np.linspace(y0, y1, n))
    q = np.vstack([X.ravel(), Y.ravel(), np.full(X.size, T)])
    Pm = sb.probes(q); return (Pm @ nod['sxx']).mean(), (Pm @ nod['syy']).mean()
p11, p12, p44 = [v * 1e-11 for v in args.pi]; pl, pt = (p11 + p12 + p44) / 2, (p11 + p12 - p44) / 2
LL, WL, gL = [v * um for v in args.resL]; LT, WT, gT = [v * um for v in args.resT]
sxL, syL = avg(a - gL - LL, a - gL, 0, WL / 2); sxT, syT = avg(a - gT - WT, a - gT, 0, LT / 2)
aL, bT = pl * sxL + pt * syL, pl * syT + pt * sxT
vout = (aL - bT) / (2 + aL + bT)
print(json.dumps(dict(model="nonlinear ideal clamp", side_um=args.side, t_um=args.t, p_kPa=args.p, nxy=args.nxy, nz=args.nz,
                      w_centre_um=w0 / um, w_over_t=w0 / T, sxL_MPa=sxL / 1e6, vout_mV_per_V=vout * 1e3,
                      newton_log=log, seconds=round(time.time() - t0, 1))))
