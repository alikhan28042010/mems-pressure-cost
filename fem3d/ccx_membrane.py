import argparse, json, os, subprocess, time
import numpy as np
import gmsh
from scipy.interpolate import LinearNDInterpolator

ap = argparse.ArgumentParser()
ap.add_argument("--side", type=float, default=600.0)
ap.add_argument("--t", type=float, default=13.0)
ap.add_argument("--p", type=float, default=100.0)
ap.add_argument("--frame", type=float, nargs=2, default=None)
ap.add_argument("--fine", type=float, default=4.0)
ap.add_argument("--coarse", type=float, default=20.0)
ap.add_argument("--pi", type=float, nargs=3, default=[6.6, -1.1, 138.1])
ap.add_argument("--res", type=float, nargs=3, default=[30, 5, 5])
ap.add_argument("--out", default="run")
args = ap.parse_args()

t0 = time.time()
os.makedirs(args.out, exist_ok=True)
a, T = args.side / 2, args.t
tol = 1e-6

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
occ = gmsh.model.occ
if args.frame:
    F, H = args.frame
    outer = occ.addBox(0, 0, -H, a + F, a + F, H + T)
    hole = occ.addBox(0, 0, -H, a, a, H)
    occ.cut([(3, outer)], [(3, hole)])
    xmax = a + F
else:
    occ.addBox(0, 0, 0, a, a, T)
    xmax = a
occ.synchronize()
vols = [v[1] for v in gmsh.model.getEntities(3)]
gmsh.model.addPhysicalGroup(3, vols, 1)

fm = gmsh.model.mesh.field
band = 60.0
b1 = fm.add("Box")
for k, v in dict(VIn=args.fine, VOut=args.coarse, XMin=a - band, XMax=a + 20, YMin=0, YMax=a + 20,
                 ZMin=-20, ZMax=T + 1, Thickness=40).items():
    fm.setNumber(b1, k, v)
b2 = fm.add("Box")
for k, v in dict(VIn=args.fine, VOut=args.coarse, XMin=0, XMax=a + 20, YMin=a - band, YMax=a + 20,
                 ZMin=-20, ZMax=T + 1, Thickness=40).items():
    fm.setNumber(b2, k, v)
b3 = fm.add("Box")
for k, v in dict(VIn=2 * args.fine, VOut=args.coarse, XMin=0, XMax=a, YMin=0, YMax=a,
                 ZMin=-1, ZMax=T + 1, Thickness=20).items():
    fm.setNumber(b3, k, v)
mn = fm.add("Min")
fm.setNumbers(mn, "FieldsList", [b1, b2, b3])
fm.setAsBackgroundMesh(mn)
for opt in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
    gmsh.option.setNumber(opt, 0)
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.option.setNumber("Mesh.HighOrderOptimize", 0)
gmsh.model.mesh.generate(3)

tags, xyz, _ = gmsh.model.mesh.getNodes()
X = xyz.reshape(-1, 3)
idx = {int(t): i for i, t in enumerate(tags)}
ets, ens = gmsh.model.mesh.getElementsByType(11)
con = np.array([idx[int(n)] for n in ens]).reshape(-1, 10)

edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
mids = X[con[:, 4:]]
order = np.zeros((len(con), 6), int)
for k, (i, j) in enumerate(edges):
    m = 0.5 * (X[con[:, i]] + X[con[:, j]])
    order[:, k] = np.argmin(np.linalg.norm(mids - m[:, None, :], axis=2), axis=1)
con_ccx = np.hstack([con[:, :4], np.take_along_axis(con[:, 4:], order, axis=1)])

top_tris = []
for dim, tag in gmsh.model.getEntities(2):
    c = gmsh.model.occ.getCenterOfMass(dim, tag)
    if abs(c[2] - T) < tol and c[0] < a + tol and c[1] < a + tol:
        tt, tn = gmsh.model.mesh.getElementsByType(9, tag)
        top_tris.append(np.array([idx[int(n)] for n in tn]).reshape(-1, 6))
gmsh.finalize()
top_tris = np.vstack(top_tris)

load = np.zeros(len(X))
P = args.p * 1e-3
p0, p1, p2 = X[top_tris[:, 0]], X[top_tris[:, 1]], X[top_tris[:, 2]]
area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
for k in (3, 4, 5):
    np.add.at(load, top_tris[:, k], P * area / 3)

x, y, z = X.T
fix1 = np.where(np.isclose(x, 0))[0]
fix2 = np.where(np.isclose(y, 0))[0]
if args.frame:
    fix_all = np.where(np.isclose(z, -args.frame[1]))[0]
    fix1 = np.union1d(fix1, np.where(np.isclose(x, xmax))[0])
    fix2 = np.union1d(fix2, np.where(np.isclose(y, xmax))[0])
else:
    fix_all = np.where(np.isclose(x, a) | np.isclose(y, a))[0]

def block(f, ids, per=8):
    ids = list(ids)
    for i in range(0, len(ids), per):
        f.write(", ".join(str(v) for v in ids[i:i + per]) + "\n")

c45 = np.cos(np.pi / 4)
inp = os.path.join(args.out, "membrane.inp")
with open(inp, "w") as f:
    f.write("*NODE, NSET=NALL\n")
    for i, (px, py, pz) in enumerate(X):
        f.write(f"{i + 1}, {px:.9g}, {py:.9g}, {pz:.9g}\n")
    f.write("*ELEMENT, TYPE=C3D10, ELSET=EALL\n")
    for e, row in enumerate(con_ccx):
        f.write(f"{e + 1}, " + ", ".join(str(v + 1) for v in row) + "\n")
    for name, ids in (("FIX1", fix1), ("FIX2", fix2), ("FIXALL", fix_all)):
        f.write(f"*NSET, NSET={name}\n")
        block(f, ids + 1)
    f.write(f"*ORIENTATION, NAME=OR1, SYSTEM=RECTANGULAR\n{c45:.12f}, {c45:.12f}, 0., {-c45:.12f}, {c45:.12f}, 0.\n")
    f.write("*MATERIAL, NAME=SI\n*ELASTIC, TYPE=ORTHO\n")
    f.write("165700., 63900., 165700., 63900., 63900., 165700., 79600., 79600.,\n79600., 0.\n")
    f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=SI, ORIENTATION=OR1\n")
    f.write("*STEP\n*STATIC\n*BOUNDARY\nFIX1, 1, 1\nFIX2, 2, 2\nFIXALL, 1, 3\n*CLOAD\n")
    for i in np.nonzero(load)[0]:
        f.write(f"{i + 1}, 3, {-load[i]:.10g}\n")
    f.write("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n")

t_mesh = time.time() - t0
run = subprocess.run(["ccx", "-i", "membrane"], cwd=args.out, capture_output=True, text=True)
if run.returncode != 0 or "Job finished" not in run.stdout:
    raise SystemExit(run.stdout[-2000:] + run.stderr[-2000:])

def read_frd(path, key, ncomp):
    out, on = {}, False
    for line in open(path):
        if line.startswith(" -4") and line.split()[1] == key:
            on = True
            continue
        if on:
            if line.startswith(" -3"):
                break
            if line.startswith(" -1"):
                nid = int(line[3:13])
                out[nid] = [float(line[13 + 12 * k:25 + 12 * k]) for k in range(ncomp)]
    ids = np.array(sorted(out))
    return ids, np.array([out[i] for i in ids])

frd = os.path.join(args.out, "membrane.frd")
nid, U = read_frd(frd, "DISP", 3)
sid, S = read_frd(frd, "STRESS", 6)
Ufull = np.zeros((len(X), 3)); Ufull[nid - 1] = U
Sfull = np.zeros((len(X), 6)); Sfull[sid - 1] = S

on_top = np.where(np.isclose(z, T))[0]
ip_x = LinearNDInterpolator(X[on_top, :2], Sfull[on_top, 0])
ip_y = LinearNDInterpolator(X[on_top, :2], Sfull[on_top, 1])
ip_w = LinearNDInterpolator(X[on_top, :2], Ufull[on_top, 2])
w = abs(float(ip_w(0, 0))) - abs(float(ip_w(a, 0)))

def avg(x0, x1, y0, y1, n=25):
    gx, gy = np.meshgrid(np.linspace(x0, x1, n), np.linspace(y0, y1, n))
    return float(np.nanmean(ip_x(gx, gy))), float(np.nanmean(ip_y(gx, gy)))

L, W, g = args.res
sxL, syL = avg(a - g - L, a - g, 0, W / 2)
sxT, syT = avg(a - g - W, a - g, 0, L / 2)
p11, p12, p44 = [v * 1e-11 for v in args.pi]
pl, pt = (p11 + p12 + p44) / 2, (p11 + p12 - p44) / 2
aL = (pl * sxL + pt * syL) * 1e6
bT = (pl * syT + pt * sxT) * 1e6
vout = (aL - bT) / (2 + aL + bT)

line_x = np.linspace(0, a - 0.01, 301)
prof = ip_x(line_x, np.zeros_like(line_x))
np.savetxt(os.path.join(args.out, "profile_top_y0.csv"), np.column_stack([line_x, prof]),
           delimiter=",", header="x_um,sigma_x_MPa", comments="", fmt="%.5g")
np.savez_compressed(os.path.join(args.out, "solution.npz"), X=X, tets=con_ccx, U=Ufull, S=Sfull, top=top_tris)

res = dict(solver="CalculiX 2.17, C3D10", mesher=f"Gmsh {gmsh.__version__}", side_um=args.side, t_um=T,
           p_kPa=args.p, frame_um=args.frame, nodes=len(X), elements=len(con), fine_um=args.fine,
           w_centre_um=float(w), sxL_MPa=sxL, syL_MPa=syL, sxT_MPa=sxT, syT_MPa=syT,
           sigma_x_2um_in_MPa=float(ip_x(a - 2, 0)), vout_mV_per_V=vout * 1e3,
           mesh_seconds=round(t_mesh, 1), total_seconds=round(time.time() - t0, 1))
json.dump(res, open(os.path.join(args.out, "result.json"), "w"), indent=1)
print(json.dumps(res))
