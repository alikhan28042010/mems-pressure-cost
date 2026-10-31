import argparse, json, time
import numpy as np
try:
    from pypardiso import spsolve
    fmt = "csr"
except ImportError:
    from scipy.sparse.linalg import spsolve
    fmt = "csc"
from skfem import (MeshHex, Basis, ElementHex2, ElementVector, BilinearForm, LinearForm,
                   FacetBasis, solve, condense)

p = argparse.ArgumentParser()
p.add_argument("--side", type=float, default=600.0)
p.add_argument("--t", type=float, default=13.0)
p.add_argument("--p", type=float, default=100.0)
p.add_argument("--nxy", type=int, default=16)
p.add_argument("--nz", type=int, default=2)
p.add_argument("--iso", type=float, nargs=2, default=None)
p.add_argument("--pi", type=float, nargs=3, default=[6.6, -1.1, 138.1])
p.add_argument("--resL", type=float, nargs=3, default=[30, 5, 5])
p.add_argument("--resT", type=float, nargs=3, default=[30, 5, 5])
p.add_argument("--depth", type=float, default=0.0)
p.add_argument("--frame", type=float, nargs=2, default=None)
p.add_argument("--map", type=str, default=None)
p.add_argument("--profile", type=str, default=None)
p.add_argument("--out", type=str, default=None)
args = p.parse_args()

t0 = time.time()
um = 1e-6
a = args.side / 2 * um
T = args.t * um
P = args.p * 1e3

def cubic_tensor(c11, c12, c44):
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if i == j and k == l:
                        C[i, j, k, l] = c11 if i == k else c12
                    elif (i == k and j == l) or (i == l and j == k):
                        if i != j:
                            C[i, j, k, l] = c44
    return C

def rotate(C, ang):
    c, s = np.cos(ang), np.sin(ang)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return np.einsum('ia,jb,kc,ld,abcd->ijkl', R, R, R, R, C)

if args.iso:
    E, nu = args.iso[0] * 1e9, args.iso[1]
    lam, mu = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
    Ct = cubic_tensor(lam + 2 * mu, lam, mu)
    material = f"isotropic E={args.iso[0]} GPa nu={nu}"
else:
    Ct = rotate(cubic_tensor(165.7e9, 63.9e9, 79.6e9), np.pi / 4)
    material = "anisotropic Si (100), edges along <110>"

s = np.linspace(0, 1, args.nxy + 1)
xs = a * np.sin(np.pi / 2 * s)
zs = np.linspace(0, T, args.nz + 1)
if args.frame:
    F, H = args.frame[0] * um, args.frame[1] * um
    nf = max(4, args.nxy // 3)
    xf = a + F * (np.linspace(0, 1, nf + 1)[1:] ** 1.6)
    xs = np.concatenate([xs, xf])
    zb = -H * (np.linspace(0, 1, nf + 1)[1:] ** 1.6)
    zs = np.concatenate([zb[::-1], zs])
    m = MeshHex.init_tensor(xs, xs, zs)
    c = m.p[:, m.t].mean(axis=1)
    cavity = np.nonzero((c[0] < a) & (c[1] < a) & (c[2] < 0))[0]
    m = m.remove_elements(cavity)
else:
    m = MeshHex.init_tensor(xs, xs, zs)
e = ElementVector(ElementHex2())
basis = Basis(m, e)

NZ_C = [(i, j, k, l, Ct[i, j, k, l]) for i in range(3) for j in range(3) for k in range(3)
        for l in range(3) if abs(Ct[i, j, k, l]) > 1e-6 * abs(Ct).max()]

@BilinearForm
def stiffness(u, v, w):
    gu, gv = u.grad, v.grad
    out = 0.0
    for i, j, k, l, c in NZ_C:
        out = out + c * gu[k, l] * gv[i, j]
    return out

@LinearForm
def load(v, w):
    return -P * v.value[2]

K = stiffness.assemble(basis)
top = m.facets_satisfying(lambda x: np.isclose(x[2], T))
f = load.assemble(FacetBasis(m, e, facets=top))

D = basis.get_dofs
if args.frame:
    xmax = xs[-1]
    outer = D(m.facets_satisfying(lambda x: np.isclose(x[2], zs[0]))).all()
    sides = np.concatenate([D(m.facets_satisfying(lambda x: np.isclose(x[0], xmax))).all('u^1'),
                            D(m.facets_satisfying(lambda x: np.isclose(x[1], xmax))).all('u^2')])
    clamp = np.concatenate([outer, sides])
else:
    clamp = D(m.facets_satisfying(lambda x: np.isclose(x[0], a) | np.isclose(x[1], a))).all()
fixed = np.unique(np.concatenate([
    D(m.facets_satisfying(lambda x: np.isclose(x[0], 0))).all('u^1'),
    D(m.facets_satisfying(lambda x: np.isclose(x[1], 0))).all('u^2'),
    clamp,
]))

def linsolve(A, b):
    return spsolve(A.asformat(fmt), b)

u = solve(*condense(K, f, D=fixed), solver=linsolve)

cdofs = basis.get_dofs(nodes=lambda x: np.isclose(x[0], 0) & np.isclose(x[1], 0)).nodal['u^3']
w_centre = float(np.mean(np.abs(u[cdofs])))
edofs = basis.get_dofs(nodes=lambda x: np.isclose(x[0], a) & np.isclose(x[1], 0) & np.isclose(x[2], T)).nodal['u^3']
w_centre -= float(np.abs(u[edofs]).mean())

sb = basis.with_element(ElementHex2())
uh = basis.interpolate(u)
strain = 0.5 * (uh.grad + np.swapaxes(uh.grad, 0, 1))
sig = np.einsum('ijkl,kl...->ij...', Ct, strain)

@BilinearForm
def mass(u_, v_, w):
    return u_ * v_

@LinearForm
def rhs(v_, w):
    return w['s'] * v_

M = mass.assemble(sb)
nodal = {}
for (i, j), name in {(0, 0): 'sxx', (1, 1): 'syy', (0, 1): 'sxy'}.items():
    b = rhs.assemble(sb, s=sig[i, j])
    nodal[name] = linsolve(M, b)

def probe(name, pts):
    return sb.probes(pts) @ nodal[name]

xline = np.linspace(0, a, 601)
pts = np.vstack([xline, np.zeros_like(xline), np.full_like(xline, T)])
prof_sxx, prof_syy = probe('sxx', pts), probe('syy', pts)
sig_edge = float(prof_sxx[-1])
sig_centre = float(prof_sxx[0])

p11, p12, p44 = [v * 1e-11 for v in args.pi]
pl, pt = (p11 + p12 + p44) / 2, (p11 + p12 - p44) / 2

def avg_rect(x0, x1, y0, y1, n=25):
    X, Y = np.meshgrid(np.linspace(x0, x1, n), np.linspace(y0, y1, n))
    zs_ = [T] if args.depth <= 0 else list(T - np.linspace(0, args.depth * um, 5))
    sx, sy = [], []
    for zz in zs_:
        q = np.vstack([X.ravel(), Y.ravel(), np.full(X.size, zz)])
        sx.append(probe('sxx', q).mean()); sy.append(probe('syy', q).mean())
    return float(np.mean(sx)), float(np.mean(sy))

LL, WL, gL = [v * um for v in args.resL]
LT, WT, gT = [v * um for v in args.resT]
sxL, syL = avg_rect(a - gL - LL, a - gL, 0, WL / 2)
dR_L = pl * sxL + pt * syL
sxT, syT = avg_rect(a - gT - WT, a - gT, 0, LT / 2)
dR_T = pl * syT + pt * sxT
vout = (dR_L - dR_T) / (2 + dR_L + dR_T)

res = dict(material=material, side_um=args.side, t_um=args.t, p_kPa=args.p, nxy=args.nxy, nz=args.nz,
           dofs=int(K.shape[0]), w_centre_um=w_centre / um, w_over_t=w_centre / T,
           sigma_edge_MPa=sig_edge / 1e6, sigma_centre_MPa=sig_centre / 1e6,
           pi_1e11=args.pi, resL_um=args.resL, resT_um=args.resT, res_depth_um=args.depth,
           frame_um=args.frame, sxL_MPa=sxL / 1e6, syL_MPa=syL / 1e6, sxT_MPa=sxT / 1e6, syT_MPa=syT / 1e6,
           sigma_x_2um_in_MPa=float(np.interp(a - 2 * um, xline, prof_sxx)) / 1e6,
           dR_L=dR_L, dR_T=dR_T, vout_mV_per_V=vout * 1e3,
           sens_mV_per_V_per_kPa=vout * 1e3 / args.p, seconds=round(time.time() - t0, 1))
print(json.dumps(res))
if args.out:
    json.dump(res, open(args.out, "w"), indent=1)
if args.map:
    g = np.linspace(0, a, 61)
    X, Y = np.meshgrid(g, g)
    q = np.vstack([X.ravel(), Y.ravel(), np.full(X.size, T)])
    np.savetxt(args.map, np.column_stack([X.ravel() / um, Y.ravel() / um, probe('sxx', q) / 1e6, probe('syy', q) / 1e6]),
               delimiter=",", header="x_um,y_um,sigma_x_MPa,sigma_y_MPa", comments="", fmt="%.5g")
if args.profile:
    np.savetxt(args.profile, np.column_stack([xline / um, prof_sxx / 1e6, prof_syy / 1e6]), delimiter=",",
               header="x_um_from_centre,sigma_x_MPa,sigma_y_MPa", comments="", fmt="%.5g")
