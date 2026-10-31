import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})

DIE, MEM = 1200, 600
c0 = DIE / 2
m0, m1 = c0 - MEM / 2, c0 + MEM / 2
PAD = 110

fig, ax = plt.subplots(figsize=(5.4, 5.4))
ax.add_patch(Rectangle((0, 0), DIE, DIE, facecolor="#c9d3df", edgecolor="#333", lw=1))
ax.add_patch(Rectangle((m0, m0), MEM, MEM, facecolor="#e8eef5", edgecolor="#555", lw=1, ls="--"))
ax.text(c0, c0, "membrane\n600 × 600 µm\n13 µm thick", ha="center", va="center", fontsize=9, color="#333")

L, Wd = 70, 18
res = {
    "top": (c0 - Wd / 2, m1 - L, Wd, L),
    "bottom": (c0 - Wd / 2, m0, Wd, L),
    "left": (m0 + 12, c0 - L / 2, Wd, L),
    "right": (m1 - 12 - Wd, c0 - L / 2, Wd, L),
}
for k, (x, y, w, h) in res.items():
    ax.add_patch(Rectangle((x, y), w, h, facecolor="#e4572e", edgecolor="#8c1c13", lw=0.8, zorder=4))

pads = {"TL": (60, DIE - 60 - PAD), "TR": (DIE - 60 - PAD, DIE - 60 - PAD), "BL": (60, 60), "BR": (DIE - 60 - PAD, 60)}
names = {"TL": "Vout+", "TR": "Vin+", "BL": "Vin−", "BR": "Vout−"}
for k, (x, y) in pads.items():
    ax.add_patch(Rectangle((x, y), PAD, PAD, facecolor="#9a9a9a", edgecolor="#444", lw=0.8, zorder=3))
    ax.text(x + PAD / 2, y + PAD / 2, names[k], ha="center", va="center", fontsize=8, color="white", zorder=5)

def wire(pts):
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color="#6d6d6d", lw=2.2, solid_capstyle="round", zorder=2)

pc = {k: (x + PAD / 2, y + PAD / 2) for k, (x, y) in pads.items()}
lx = m0 + 12 + Wd / 2
rx = m1 - 12 - Wd / 2
lo, hi = DIE - 60 - PAD, 60 + PAD
wire([(c0, m1), (c0 - 30, m1), (c0 - 30, pc["TL"][1]), (hi, pc["TL"][1])])
wire([(lx, c0 + L / 2), (lx, m1 + 60), (pc["TL"][0], m1 + 60), (pc["TL"][0], lo)])
wire([(c0, m1 - L), (c0 + 30, m1 - L), (c0 + 30, pc["TR"][1]), (lo, pc["TR"][1])])
wire([(rx, c0 + L / 2), (rx, m1 + 60), (pc["TR"][0], m1 + 60), (pc["TR"][0], lo)])
wire([(c0, m0), (c0 - 30, m0), (c0 - 30, pc["BL"][1]), (hi, pc["BL"][1])])
wire([(lx, c0 - L / 2), (lx, m0 - 60), (pc["BL"][0], m0 - 60), (pc["BL"][0], hi)])
wire([(c0, m0 + L), (c0 + 30, m0 + L), (c0 + 30, pc["BR"][1]), (lo, pc["BR"][1])])
wire([(rx, c0 - L / 2), (rx, m0 - 60), (pc["BR"][0], m0 - 60), (pc["BR"][0], hi)])

ax.annotate("resistors across the edge\n(top, bottom)", xy=(c0 - 12, m1 - L / 2), xytext=(c0 - 250, m1 - 170),
            fontsize=8, arrowprops=dict(arrowstyle="-", color="#333", lw=0.7))
ax.annotate("resistors along the edge\n(left, right)", xy=(rx - 10, c0), xytext=(m1 - 430, c0 - 150),
            fontsize=8, arrowprops=dict(arrowstyle="-", color="#333", lw=0.7))
ax.plot([0, DIE], [-45, -45], color="#333", lw=0.8)
ax.text(c0, -95, "1.2 mm", ha="center", fontsize=9)
ax.set_xlim(-30, DIE + 30)
ax.set_ylim(-130, DIE + 30)
ax.set_aspect("equal")
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig11_chip_layout.png"), dpi=250)
plt.close(fig)

steps = [
    ("Front end\n6 masks, SON", "$0.095 per die"),
    ("Electrical\nwafer test", "bridge resistance"),
    ("Dicing\n1.2 × 1.2 mm", "9,390 dies/wafer"),
    ("Readout chip\n(ASIC)", "$0.15"),
    ("Moulded\npackage", "$0.10"),
    ("Calibration\n3 p × 3 T", "$0.23"),
    ("Final test", "$0.60 per sensor"),
]
fig, ax = plt.subplots(figsize=(10, 2.2))
bw, gap = 1.25, 0.32
for i, (a, b) in enumerate(steps):
    x = i * (bw + gap)
    ax.add_patch(FancyBboxPatch((x, 0.35), bw, 1.0, boxstyle="round,pad=0.02,rounding_size=0.08",
                                facecolor="#e8eef5" if i else "#b4c3d6", edgecolor="#4a5d73", lw=1))
    ax.text(x + bw / 2, 0.85, a, ha="center", va="center", fontsize=9)
    ax.text(x + bw / 2, 0.1, b, ha="center", va="center", fontsize=8.5, color="#333")
    if i < len(steps) - 1:
        ax.annotate("", xy=(x + bw + gap - 0.03, 0.85), xytext=(x + bw + 0.03, 0.85),
                    arrowprops=dict(arrowstyle="-|>", color="#4a5d73", lw=1.2))
ax.set_xlim(-0.1, len(steps) * (bw + gap))
ax.set_ylim(-0.1, 1.45)
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig12_backend_flow.png"), dpi=250)
plt.close(fig)
print("figures/fig11_chip_layout.png figures/fig12_backend_flow.png")
