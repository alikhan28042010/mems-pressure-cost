import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

HERE = os.path.dirname(os.path.abspath(__file__))
SI, EPI, OX, NIT, PRES, PPLUS, AL, VOID = "#c9d3df", "#b4c3d6", "#f3df8a", "#7fbf7f", "#e4572e", "#8c1c13", "#6d6d6d", "white"
W, BASE, EPI_TOP = 10.0, 3.0, 4.3
CAV = (3.0, 7.0, 2.25, 2.7)
RES = [(3.05, 3.55), (6.45, 6.95)]
CON = [(2.55, 3.05), (6.95, 7.45)]

def rect(ax, x0, x1, y0, y1, c, z=1, lw=0.0):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=c, edgecolor="#333" if lw else c, lw=lw, zorder=z))

def wafer(ax, top):
    rect(ax, 0, W, 0, BASE, SI)
    if top > BASE:
        rect(ax, 0, W, BASE, top, EPI)

def cavity(ax):
    rect(ax, CAV[0], CAV[1], CAV[2], CAV[3], VOID, 3)

def trenches(ax):
    for k in range(9):
        x = 3.15 + k * 0.45
        rect(ax, x, x + 0.2, 2.0, BASE, VOID, 3)

def oxide(ax, top, t=0.12, holes=()):
    rect(ax, 0, W, top, top + t, OX, 4)
    for x0, x1 in holes:
        rect(ax, x0, x1, top, top + t, VOID, 5)

def resistors(ax, top):
    for x0, x1 in RES:
        rect(ax, x0, x1, top - 0.18, top, PRES, 5)

def contacts(ax, top):
    for x0, x1 in CON:
        rect(ax, x0, x1, top - 0.3, top, PPLUS, 5)

steps = [
    ("a  Oxide mask on n-type (100) wafer", lambda ax: (wafer(ax, BASE), oxide(ax, BASE, holes=[(3.0, 7.0)]))),
    ("b  DRIE trench array, 0.8 µm holes", lambda ax: (wafer(ax, BASE), trenches(ax))),
    ("c  H₂ anneal 1100 °C, 15 min: empty space\n    forms under ~1 µm of silicon", lambda ax: (wafer(ax, BASE), cavity(ax))),
    ("d  Epitaxy 12 µm: membrane 13 µm", lambda ax: (wafer(ax, EPI_TOP), cavity(ax))),
    ("e  Thin oxide, boron implant:\n    piezoresistors 200–500 Ω/sq", lambda ax: (wafer(ax, EPI_TOP), cavity(ax), resistors(ax, EPI_TOP), oxide(ax, EPI_TOP, 0.06))),
    ("f  P+ implant for contacts", lambda ax: (wafer(ax, EPI_TOP), cavity(ax), resistors(ax, EPI_TOP), contacts(ax, EPI_TOP), oxide(ax, EPI_TOP, 0.06))),
    ("g  SiO₂ + Si₃N₄, contact holes", lambda ax: (wafer(ax, EPI_TOP), cavity(ax), resistors(ax, EPI_TOP), contacts(ax, EPI_TOP),
                                                     oxide(ax, EPI_TOP, 0.14, holes=[(2.65, 2.95), (7.05, 7.35)]),
                                                     rect(ax, 0, 2.65, EPI_TOP + 0.14, EPI_TOP + 0.26, NIT, 4),
                                                     rect(ax, 2.95, 7.05, EPI_TOP + 0.14, EPI_TOP + 0.26, NIT, 4),
                                                     rect(ax, 7.35, W, EPI_TOP + 0.14, EPI_TOP + 0.26, NIT, 4))),
    ("h  Al interconnect: Wheatstone bridge", lambda ax: (wafer(ax, EPI_TOP), cavity(ax), resistors(ax, EPI_TOP), contacts(ax, EPI_TOP),
                                                          oxide(ax, EPI_TOP, 0.14, holes=[(2.65, 2.95), (7.05, 7.35)]),
                                                          rect(ax, 0, 2.65, EPI_TOP + 0.14, EPI_TOP + 0.26, NIT, 4),
                                                          rect(ax, 2.95, 7.05, EPI_TOP + 0.14, EPI_TOP + 0.26, NIT, 4),
                                                          rect(ax, 7.35, W, EPI_TOP + 0.14, EPI_TOP + 0.26, NIT, 4),
                                                          rect(ax, 1.2, 3.0, EPI_TOP + 0.14, EPI_TOP + 0.45, AL, 6),
                                                          rect(ax, 2.65, 2.95, EPI_TOP, EPI_TOP + 0.14, AL, 6),
                                                          rect(ax, 7.0, 8.8, EPI_TOP + 0.14, EPI_TOP + 0.45, AL, 6),
                                                          rect(ax, 7.05, 7.35, EPI_TOP, EPI_TOP + 0.14, AL, 6))),
]

fig, axes = plt.subplots(2, 4, figsize=(13, 5.4))
for ax, (label, draw) in zip(axes.ravel(), steps):
    draw(ax)
    ax.set_xlim(0, W)
    ax.set_ylim(0, 5.0)
    ax.set_aspect("auto")
    ax.axis("off")
    ax.set_title(label, fontsize=9.5, loc="left")
handles = [Patch(color=c, label=n) for c, n in ((SI, "n-type Si substrate"), (EPI, "epitaxial Si"), (OX, "SiO₂"),
                                                (NIT, "Si₃N₄"), (PRES, "p piezoresistor"), (PPLUS, "p+ contact"),
                                                (AL, "aluminium"))]
handles.append(Patch(facecolor="white", edgecolor="#333", label="vacuum cavity / trench"))
fig.legend(handles=handles, loc="lower center", ncol=8, frameon=False, fontsize=9)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(os.path.join(HERE, "figures", "fig10_process_flow.png"), dpi=200)
print("figures/fig10_process_flow.png")
