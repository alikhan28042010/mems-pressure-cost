import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
res3d = os.path.join(ROOT, "results", "fem3d")
figs = os.path.join(ROOT, "figures")

ccx = list(csv.DictReader(open(os.path.join(res3d, "ccx_runs.csv"))))
sk_frame = {float(r["t_um"]): r for r in csv.DictReader(open(os.path.join(ROOT, "data", "sweep_su.csv")))}
sk_clamp = {float(r["t_um"]): r for r in csv.DictReader(open(os.path.join(ROOT, "data", "sweep_su_clamp.csv")))}

best = {}
for r in ccx:
    k = (float(r["t_um"]), r["frame_um"])
    if k not in best or float(r["fine_um"]) < float(best[k]["fine_um"]):
        best[k] = r

rows = []
for (t, fr), r in sorted(best.items()):
    sk = (sk_frame if fr != "none" else sk_clamp)[t]
    w1, w2 = float(sk["w_centre_um"]), float(r["w_centre_um"])
    v1, v2 = float(sk["vout_mV_per_V"]), float(r["vout_mV_per_V"])
    rows.append(dict(t_um=t, edge="frame" if fr != "none" else "clamp", ccx_fine_um=r["fine_um"], ccx_nodes=r["nodes"],
                     w_skfem_um=round(w1, 4), w_ccx_um=round(w2, 4), w_diff_pct=round(100 * (w2 / w1 - 1), 2),
                     vout_skfem=round(v1, 3), vout_ccx=round(v2, 3), vout_diff_pct=round(100 * (v2 / v1 - 1), 2)))
with open(os.path.join(res3d, "crosscheck_skfem_vs_calculix.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)

fr = [r for r in rows if r["edge"] == "frame"]
t = np.array([r["t_um"] for r in fr])
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
ts = sorted(sk_frame)
ax[0].plot(ts, [float(sk_frame[k]["vout_mV_per_V"]) for k in ts], "-", color="#2a6fdb", label="scikit-fem, hex27")
ax[0].plot(t, [r["vout_ccx"] for r in fr], "o", mfc="none", color="#d9480f", ms=7, label="CalculiX, tet10")
ax[0].set_xlabel("Membrane thickness (µm)"); ax[0].set_ylabel("Bridge output at 100 kPa (mV/V)")
ax[1].plot(ts, [float(sk_frame[k]["w_centre_um"]) for k in ts], "-", color="#2a6fdb", label="scikit-fem, hex27")
ax[1].plot(t, [r["w_ccx_um"] for r in fr], "o", mfc="none", color="#d9480f", ms=7, label="CalculiX, tet10")
ax[1].set_xlabel("Membrane thickness (µm)"); ax[1].set_ylabel("Centre deflection (µm)")
for a in ax:
    a.grid(alpha=0.3); a.legend(frameon=False, fontsize=9)
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(figs, "fig8_solver_check.png"), dpi=200)
plt.close(fig)

d = np.load(os.path.join(res3d, "solution_13um_frame.npz"))
X, U, S, top = d["X"], d["U"], d["S"], d["top"]
tri = top[:, :3]
a_half = 300.0
keep = np.all(X[tri][:, :, 0] <= a_half + 1e-6, axis=1) & np.all(X[tri][:, :, 1] <= a_half + 1e-6, axis=1)
tri = tri[keep]
scale = 150.0
fig = plt.figure(figsize=(7, 5.2))
ax = fig.add_subplot(projection="3d")
norm = colors.TwoSlopeNorm(vmin=-30, vcenter=0, vmax=60)
for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
    P = X + scale * U
    xx, yy, zz = sx * P[:, 0], sy * P[:, 1], P[:, 2]
    val = S[tri, 0].mean(axis=1)
    poly = ax.plot_trisurf(xx, yy, zz, triangles=tri, linewidth=0, antialiased=False, shade=False)
    poly.set_facecolor(cm.RdBu_r(norm(val)))
ax.set_box_aspect((1, 1, 0.35))
ax.view_init(elev=28, azim=-55)
ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)"); ax.set_zlabel(f"z (µm, ×{scale:.0f})")
sm = cm.ScalarMappable(norm=norm, cmap="RdBu_r")
fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.16, label="σx on top surface (MPa)")
fig.tight_layout()
fig.savefig(os.path.join(figs, "fig9_membrane_3d.png"), dpi=200)
print("figures saved")
