import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, FIG = os.path.join(HERE, "data"), os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
BLUE, ORANGE, AQUA, VIOLET, RED, GOLD = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#d6455d", "#b8860b"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.axisbelow": True, "figure.facecolor": "white", "savefig.dpi": 200})

def read(name):
    return list(csv.DictReader(open(os.path.join(DATA, name))))

fig, ax = plt.subplots(figsize=(7.2, 3.8))
cmap = plt.get_cmap("viridis")
ts = [9, 11, 13, 15, 17]
for i, t in enumerate(ts):
    f = os.path.join(DATA, f"profile_frame_t{t}.csv")
    if not os.path.exists(f): continue
    r = np.loadtxt(f, delimiter=",", skiprows=1)
    ax.plot(300 - r[:, 0], r[:, 1], color=cmap(i / (len(ts) - 1) * 0.9), lw=2, label=f"{t} µm")
ax.axvspan(5, 35, color=GOLD, alpha=0.12, lw=0)
ax.text(20, ax.get_ylim()[1] * 0.92, "resistor", ha="center", color=MUTED, fontsize=9)
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xlim(0, 300); ax.set_xlabel("Distance from the membrane edge (µm)")
ax.set_ylabel("Top-surface stress σ$_x$ (MPa)")
ax.legend(title="Thickness", frameon=False, ncol=1, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_stress_profile.png")); plt.close(fig)

fr, cl = read("sweep_su.csv"), read("sweep_su_clamp.csv")
fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.6))
for rows, col, lab, mk in ((fr, BLUE, "frame model", "o"), (cl, ORANGE, "ideal clamp", "s")):
    t = np.array([float(r["t_um"]) for r in rows])
    v = np.array([float(r["vout_mV_per_V"]) for r in rows]); w = np.array([float(r["w_centre_um"]) for r in rows])
    nv = -np.polyfit(np.log(t), np.log(v), 1)[0]; nw = -np.polyfit(np.log(t), np.log(w), 1)[0]
    ax[0].loglog(t, v, mk + "-", color=col, ms=5, lw=1.5, label=f"{lab} (∝ t$^{{-{nv:.2f}}}$)")
    ax[1].loglog(t, w, mk + "-", color=col, ms=5, lw=1.5, label=f"{lab} (∝ t$^{{-{nw:.2f}}}$)")
ax[1].annotate("Su et al. (2019), ANSYS: 0.64 µm", xy=(13, 0.64), xytext=(13.4, 0.9), color=MUTED, fontsize=8,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax[1].plot([13], [0.64], "D", color=MUTED, ms=5)
for a, yl in zip(ax, ("Bridge output at 100 kPa (mV/V)", "Centre deflection at 100 kPa (µm)")):
    a.set_xticks([9, 11, 13, 15, 17]); a.set_xticklabels(["9", "11", "13", "15", "17"]); a.minorticks_off()
    a.set_xlabel("Membrane thickness t (µm)"); a.set_ylabel(yl); a.legend(frameon=False, fontsize=8)
ax[0].set_yticks([20, 30, 50, 70]); ax[0].set_yticklabels(["20", "30", "50", "70"])
ax[1].set_yticks([0.3, 0.5, 1, 1.5]); ax[1].set_yticklabels(["0.3", "0.5", "1", "1.5"])
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_thickness_law.png")); plt.close(fig)

sc = read("singh_comparison.csv")
def pick(edge, nxy, depth):
    return sorted([(float(r["t_um"]), float(r["this_model_mV_per_V_per_bar"])) for r in sc
                   if r["edge_condition"] == edge and r["nxy"] == nxy and abs(float(r["res_depth_um"]) - depth) < 1e-9])
fig, ax = plt.subplots(figsize=(7.2, 3.8))
for r in sc:
    if r["singh_measured_low"] and r["edge_condition"] == "frame" and r["nxy"] == "12" and float(r["res_depth_um"]) == 1.0:
        t = float(r["t_um"])
        ax.add_patch(plt.Rectangle((t - 1, float(r["singh_measured_low"])), 2, float(r["singh_measured_high"]) - float(r["singh_measured_low"]),
                                   color=AQUA, alpha=0.25, lw=0))
ax.plot([], [], "s", color=AQUA, alpha=0.4, ms=9, label="Singh et al. (2014), measured (t \u00b1 1 \u00b5m)")
sim = sorted({(float(r["t_um"]), float(r["singh_simulated"])) for r in sc if r["singh_simulated"]})
ax.plot([s_[0] for s_ in sim], [s_[1] for s_ in sim], "s--", color=ORANGE, ms=5, lw=1.3, label="Singh et al. (2014), FEM")
m0 = pick("frame", "12", 0.0); m1 = pick("frame", "12", 1.0)
ax.plot([m[0] for m in m0], [m[1] for m in m0], "o:", color=MUTED, ms=4, lw=1.2, label="This model, surface stress only")
ax.plot([m[0] for m in m1], [m[1] for m in m1], "o-", color=BLUE, ms=6, lw=2, label="This model, averaged over 1 \u00b5m resistor depth")
ax.set_xlabel("Membrane thickness (\u00b5m)"); ax.set_ylabel("Sensitivity (mV/V/bar)")
ax.legend(frameon=False, fontsize=8.5); fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig4_singh_validation.png")); plt.close(fig)

rows = read("process_curve_base.csv"); names = [k for k in rows[0] if k != "G_pct" and not k.startswith("Silicon-on-nothing")]
LBL = {n: n for n in names}; LBL["Cavity-SOI, EC-SOI (±0.2 µm)"] = "EC-SOI or silicon-on-nothing (±0.2 µm)"
cols = [BLUE, ORANGE, AQUA, GOLD, VIOLET, RED]
fig, ax = plt.subplots(figsize=(7.4, 3.9))
G = np.array([float(r["G_pct"]) for r in rows])
styles = ["-", "-", "-", "--", "-", ":"]
for n, c, ls in zip(names, cols, styles):
    ax.plot(G, [100 * float(r[n]) for r in rows], ls, color=c, lw=2, label=LBL[n])
ax.set_xscale("log"); ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(["±1", "±2", "±5", "±10", "±20"]); ax.minorticks_off()
ax.set_xlabel("Gain range the readout chip can correct (% of nominal sensitivity)"); ax.set_ylabel("Chips inside the range (%)")
ax.set_ylim(0, 102); ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig5_process_yield.png")); plt.close(fig)

cb = read("cost_breakeven.csv")
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.7), sharey=True)
for a, screen, title in zip(ax, ("final", "probe"), ("Bad chips found at final test", "Bad chips removed at wafer probe")):
    for n, c, ls in zip(names, cols, styles):
        if n.startswith("Timed"): continue
        pts = sorted([(float(r["G_pct"]), float(r["breakeven_wafer_multiple"])) for r in cb if r["screen"] == screen and r["route"] == n])
        a.plot([p[0] for p in pts], [p[1] for p in pts], "o" + ls, color=c, lw=1.8, ms=4, label=LBL[n])
    a.axhline(1, color=MUTED, lw=1, ls="--"); a.set_yscale("log"); a.set_title(title, fontsize=10, loc="left")
    a.set_xticks([5, 10, 15, 20, 25]); a.set_xticklabels(["±5", "±10", "±15", "±20", "±25"])
    a.set_xlabel("Correctable gain range (%)"); a.minorticks_off()
ax[0].set_yticks([1, 1.5, 2, 3, 5, 10]); ax[0].set_yticklabels(["1×", "1.5×", "2×", "3×", "5×", "10×"]); ax[0].set_ylim(0.9, 12)
ax[0].set_ylabel("Break-even wafer price\n(× timed-etch wafer)"); ax[1].legend(frameon=False, fontsize=7.5)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig6_breakeven.png")); plt.close(fig)

tr = [r for r in read("cost_tornado.csv") if r["route"].startswith("Timed")]
labels = {"C_wafer": "Wafer price", "Y_wafer": "Wafer yield", "C_asic": "Readout chip price", "C_pkg": "Package price",
          "rate_cal": "Calibration chamber cost", "t_cal": "Calibration time", "Y_other": "Assembly yield",
          "C_probe": "Wafer-probe test price", "G": "Correctable gain range"}
tr = [r for r in tr if r["input"] != "C_probe"]
tr.sort(key=lambda r: abs(float(r["cost_at_high"]) - float(r["cost_at_low"])))
base = float(tr[0]["base_cost"])
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for i, r in enumerate(tr):
    lo, hi = float(r["cost_at_low"]), float(r["cost_at_high"])
    ax.barh(i, lo - base, left=base, color=BLUE, height=0.6)
    ax.barh(i, hi - base, left=base, color=ORANGE, height=0.6)
ax.set_yticks(range(len(tr))); ax.set_yticklabels([labels[r["input"]] for r in tr])
ax.axvline(base, color=INK, lw=1); ax.set_xlabel("Cost per calibrated sensor, timed wet etch (US$)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=BLUE, label="input at its low value"), Patch(color=ORANGE, label="input at its high value")], frameon=False, fontsize=8, loc="lower right"); ax.grid(axis="y", visible=False)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig7_tornado.png")); plt.close(fig)

f = os.path.join(DATA, "map_frame_t13.csv")
if os.path.exists(f):
    m = np.loadtxt(f, delimiter=",", skiprows=1)
    n = int(round(np.sqrt(len(m))))
    X, Y, S = m[:, 0].reshape(n, n), m[:, 1].reshape(n, n), m[:, 2].reshape(n, n)
    fig, ax = plt.subplots(figsize=(4.6, 3.9))
    lim = np.abs(S).max()
    cs = ax.contourf(X, Y, S, levels=np.linspace(-lim, lim, 25), cmap="RdBu_r")
    cb_ = fig.colorbar(cs, ax=ax); cb_.set_label("σ$_x$ (MPa)")
    ax.add_patch(plt.Rectangle((300 - 35, 0), 30, 2.5, fill=False, ec=INK, lw=1))
    ax.set_xlabel("x from centre (µm)"); ax.set_ylabel("y from centre (µm)"); ax.set_aspect("equal"); ax.grid(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_stress_map.png")); plt.close(fig)
print("figures written to", FIG)
