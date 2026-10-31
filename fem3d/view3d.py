import argparse, json
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("--scale", type=float, default=150.0)
ap.add_argument("--html", default="membrane_3d.html")
ap.add_argument("--vtu", default="membrane_3d.vtu")
args = ap.parse_args()

d = np.load(f"{args.run}/solution.npz")
X, tets, U, S = d["X"], d["tets"], d["U"], d["S"]
info = json.load(open(f"{args.run}/result.json"))

faces = np.vstack([tets[:, [0, 2, 1]], tets[:, [0, 1, 3]], tets[:, [1, 2, 3]], tets[:, [0, 3, 2]]])
key = np.sort(faces, axis=1)
_, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
skin = faces[cnt[inv.ravel()] == 1]
sym = np.all(np.isclose(X[skin][:, :, 0], 0), axis=1) | np.all(np.isclose(X[skin][:, :, 1], 0), axis=1)
skin = skin[~sym]

sx, sy, sz, sxy, sxz, syz = S.T
mises = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3 * (sxy ** 2 + sxz ** 2 + syz ** 2))

pts, tris, val = [], [], []
n = len(X)
for k, (fx, fy) in enumerate(((1, 1), (-1, 1), (-1, -1), (1, -1))):
    P = X + U * args.scale
    P = P * np.array([fx, fy, 1.0])
    pts.append(P)
    t = skin + k * n
    if fx * fy < 0:
        t = t[:, [0, 2, 1]]
    tris.append(t)
    val.append(sx)
pts, tris, val = np.vstack(pts), np.vstack(tris), np.concatenate(val)

used = np.unique(tris)
remap = -np.ones(len(pts), int)
remap[used] = np.arange(len(used))
pts, tris, val = pts[used], remap[tris], val[used]

frame = info["frame_um"]
title = (f"{info['t_um']:.0f} µm membrane, {info['p_kPa']:.0f} kPa, "
         f"{'membrane on frame' if frame else 'ideal clamp'} | CalculiX {info['elements']} C3D10 elements | "
         f"centre deflection {info['w_centre_um']:.3f} µm (shown ×{args.scale:.0f})")

import plotly.graph_objects as go
fig = go.Figure(go.Mesh3d(x=pts[:, 0], y=pts[:, 1], z=pts[:, 2], i=tris[:, 0], j=tris[:, 1], k=tris[:, 2],
                          intensity=val, colorscale="RdBu_r", cmid=0, cmin=-60, cmax=60,
                          colorbar=dict(title="σx (MPa)"), flatshading=False,
                          lighting=dict(ambient=0.55, diffuse=0.7, specular=0.15),
                          hovertemplate="x %{x:.0f} µm<br>y %{y:.0f} µm<br>σx %{intensity:.1f} MPa<extra></extra>"))
fig.update_layout(title=dict(text=title, font=dict(size=13)), margin=dict(l=0, r=0, t=40, b=0),
                  scene=dict(aspectmode="manual", aspectratio=dict(x=1, y=1, z=0.35),
                             xaxis_title="x (µm)", yaxis_title="y (µm)", zaxis_title="z (µm, deflection ×%d)" % args.scale,
                             camera=dict(eye=dict(x=1.4, y=-1.6, z=1.1))))
fig.write_html(args.html, include_plotlyjs=False, full_html=True)
html = open(args.html).read().replace(
    "<head>", '<head><meta charset="utf-8"><title>MEMS membrane 3D</title>'
    '<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>', 1)
open(args.html, "w").write(html)

import meshio
cells = [("tetra10", tets)]
meshio.write(args.vtu, meshio.Mesh(X, cells, point_data={"U": U, "sigma_x": sx, "sigma_y": sy, "von_mises": mises}))
print(args.html, len(pts), "points", len(tris), "triangles")
