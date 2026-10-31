import numpy as np, csv, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)
rng = np.random.default_rng(2026)
N_MC = 200_000
T0 = 13.0

rows = list(csv.DictReader(open(os.path.join(DATA, "sweep_su.csv"))))
t_tab = np.array([float(r["t_um"]) for r in rows])
v_tab = np.array([float(r["vout_mV_per_V"]) for r in rows])
order = np.argsort(t_tab); t_tab, v_tab = t_tab[order], v_tab[order]
logfit = np.polyfit(np.log(t_tab), np.log(v_tab), 1)

def sens(t):
    t = np.asarray(t, float)
    inside = (t >= t_tab[0]) & (t <= t_tab[-1])
    out = np.exp(np.polyval(logfit, np.log(np.clip(t, 1e-3, None))))
    out[inside] = np.exp(np.interp(np.log(t[inside]), np.log(t_tab), np.log(v_tab)))
    return out

S0 = float(sens(np.array([T0]))[0])

PUB = os.path.join(DATA, "public")
ROUTES, TOL, DIE = {}, {}, {}
for r in csv.DictReader(open(os.path.join(PUB, "thickness_control.csv"), encoding="utf-8")):
    ROUTES[r["route"]] = (float(r["sigma_base_um"]), float(r["sigma_conservative_um"]), r["source"])
    TOL[r["route"]] = float(r["band_um"])
    DIE[r["route"]] = float(r["die_mm"])
UNPUBLISHED = {}
for r in csv.DictReader(open(os.path.join(PUB, "thickness_scenarios.csv"), encoding="utf-8")):
    UNPUBLISHED.setdefault(r["route"], []).append(float(r["sigma_um"]))
GAINS = [2, 5, 10, 15, 20, 25]

def share_in(sig, G):
    t = T0 + rng.normal(0, sig, N_MC)
    s = sens(t) / S0
    return float(np.mean(np.abs(s - 1) <= G / 100)), float(np.std(s) / np.mean(s))

proc_rows = []
for name, (sb, sc, src) in ROUTES.items():
    for case, sig in (("base", sb), ("conservative", sc)):
        r = {"route": name, "case": case, "sigma_um": round(sig, 4), "source": src}
        for G in GAINS:
            sh, cv = share_in(sig, G)
            r[f"in_pm{G}pct"] = round(sh, 4)
        r["sensitivity_CV_pct"] = round(100 * cv, 2)
        proc_rows.append(r)
for name, sigs in UNPUBLISHED.items():
    for sig in sigs:
        r = {"route": name, "case": f"sigma={sig}", "sigma_um": sig, "source": "ASSUMPTION"}
        for G in GAINS:
            sh, cv = share_in(sig, G)
            r[f"in_pm{G}pct"] = round(sh, 4)
        r["sensitivity_CV_pct"] = round(100 * cv, 2)
        proc_rows.append(r)

with open(os.path.join(DATA, "process_montecarlo.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(proc_rows[0].keys())); w.writeheader(); w.writerows(proc_rows)

Gs = np.round(np.geomspace(0.5, 30, 40), 3)
curve = [["G_pct"] + [n for n in ROUTES]]
samples = {n: sens(T0 + np.random.default_rng(7).normal(0, ROUTES[n][0], N_MC)) / S0 for n in ROUTES}
for G in Gs:
    curve.append([G] + [round(float(np.mean(np.abs(samples[n] - 1) <= G / 100)), 4) for n in ROUTES])
with open(os.path.join(DATA, "process_curve_base.csv"), "w", newline="") as fh:
    csv.writer(fh).writerows(curve)

INPUTS = {}
for r in csv.DictReader(open(os.path.join(PUB, "cost_model_inputs.csv"), encoding="utf-8")):
    INPUTS[r["name"]] = (float(r["base"]), float(r["low"]), float(r["high"]), r["source"] or "ASSUMPTION")
def n_gross(die_mm, street=0.1, wafer=150.0, edge=3.0):
    D, A = wafer - 2 * edge, (die_mm + street) ** 2
    return int(np.pi * (D / 2) ** 2 / A - np.pi * D / np.sqrt(2 * A))
N_GROSS = n_gross(1.2)

def y_route(sig, G):
    t = T0 + np.random.default_rng(11).normal(0, sig, N_MC)
    return float(np.mean(np.abs(sens(t) / S0 - 1) <= G / 100))

def cost(p, sig, screen="final", wafer_mult=1.0, N=None):
    yg = max(y_route(sig, p["G"]), 1e-6)
    c_cal = p["rate_cal"] * p["t_cal"] / 3600
    silicon = p["C_wafer"] * wafer_mult / ((N or p["N"]) * p["Y_wafer"])
    if screen == "final":
        c = (silicon + p["C_asic"] + p["C_pkg"] + c_cal) / yg
    else:
        c = (silicon + p["C_probe"]) / yg + p["C_asic"] + p["C_pkg"] + c_cal
    return c / p["Y_other"]

base = {k: v[0] for k, v in INPUTS.items()}; base["N"] = N_GROSS
REF = "Timed wet etch (KOH)"

def breakeven(route, p, screen, die=False):
    sig_ref, sig = ROUTES[REF][0], ROUTES[route][0]
    n_ref = n_gross(DIE[REF]) if die else None
    n_r = n_gross(DIE[route]) if die else None
    c_ref = cost(p, sig_ref, screen, N=n_ref)
    lo, hi = 0.1, 100.0
    if cost(p, sig, screen, hi, N=n_r) <= c_ref: return hi
    for _ in range(60):
        mid = (lo + hi) / 2
        if cost(p, sig, screen, mid, N=n_r) > c_ref: hi = mid
        else: lo = mid
    return lo

cost_rows = []
for G in [5, 10, 15, 20, 25]:
    p = dict(base); p["G"] = G
    for screen in ("final", "probe"):
        for route in ROUTES:
            sig = ROUTES[route][0]
            cost_rows.append({"G_pct": G, "screen": screen, "route": route,
                              "yield_in_range": round(y_route(sig, G), 4),
                              "cost_same_wafer_price_usd": round(cost(p, sig, screen), 4),
                              "breakeven_wafer_multiple": round(breakeven(route, p, screen), 3) if route != REF else 1.0,
                              "die_mm": DIE[route], "dies_per_wafer": n_gross(DIE[route]),
                              "breakeven_with_route_die_size": round(breakeven(route, p, screen, die=True), 3) if route != REF else 1.0})
with open(os.path.join(DATA, "cost_breakeven.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(cost_rows[0].keys())); w.writeheader(); w.writerows(cost_rows)

p = dict(base); yg = y_route(ROUTES[REF][0], p["G"])
silicon = p["C_wafer"] / (p["N"] * p["Y_wafer"])
back = p["C_asic"] + p["C_pkg"] + p["rate_cal"] * p["t_cal"] / 3600
summary = {"N_gross": N_GROSS, "N_gross_2mm_die": n_gross(2.0), "silicon_cost_per_die_usd": round(silicon, 4),
           "backend_cost_per_die_usd": round(back, 4), "backend_share": round(back / (silicon + back), 3),
           "calibration_cost_usd": round(p["rate_cal"] * p["t_cal"] / 3600, 4), "S0_mV_per_V": round(S0, 3),
           "fem_sweep_logfit_exponent": round(-logfit[0], 3)}

torn = []
for route in (REF, "Cavity-SOI, standard (±0.5 µm)"):
    c0 = cost(base, ROUTES[route][0])
    for k, (b, lo, hi, src) in INPUTS.items():
        vals = []
        for v in (lo, hi):
            p = dict(base); p[k] = v
            vals.append(cost(p, ROUTES[route][0]))
        torn.append({"route": route, "input": k, "base_cost": round(c0, 4), "cost_at_low": round(vals[0], 4),
                     "cost_at_high": round(vals[1], 4), "low": lo, "high": hi, "source": src})
with open(os.path.join(DATA, "cost_tornado.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(torn[0].keys())); w.writeheader(); w.writerows(torn)
with open(os.path.join(DATA, "cost_inputs.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["input", "base", "low", "high", "source"])
    for k, v in INPUTS.items(): w.writerow([k, *v])
def y_general(route, G, dop_cv=0.0, shape="normal", seed=11):
    r = np.random.default_rng(seed)
    if shape == "normal": t = T0 + r.normal(0, ROUTES[route][0], N_MC)
    else: t = T0 + r.uniform(-TOL[route], TOL[route], N_MC)
    s = sens(t) / S0 * (1 + r.normal(0, dop_cv, N_MC))
    return float(np.mean(np.abs(s - 1) <= G / 100))

def cost_y(p, yg, screen, wafer_mult=1.0, N=None):
    c_cal = p["rate_cal"] * p["t_cal"] / 3600
    silicon = p["C_wafer"] * wafer_mult / ((N or p["N"]) * p["Y_wafer"])
    if screen == "final": c = (silicon + p["C_asic"] + p["C_pkg"] + c_cal) / yg
    else: c = (silicon + p["C_probe"]) / yg + p["C_asic"] + p["C_pkg"] + c_cal
    return c / p["Y_other"]

def breakeven_y(y_ref, y_r, p, screen):
    c_ref = cost_y(p, y_ref, screen); lo, hi = 0.1, 100.0
    if cost_y(p, y_r, screen, hi) <= c_ref: return hi
    for _ in range(60):
        mid = (lo + hi) / 2
        if cost_y(p, y_r, screen, mid) > c_ref: hi = mid
        else: lo = mid
    return lo

rob = []
for label, kw in (("normal, no doping spread", {}), ("normal + doping CV 2 %", {"dop_cv": 0.02}),
                  ("normal + doping CV 4 %", {"dop_cv": 0.04}), ("uniform over full tolerance", {"shape": "uniform"})):
    for G in (5, 10, 20):
        p = dict(base); p["G"] = G
        yref = y_general(REF, G, **kw)
        for route in ROUTES:
            y = y_general(route, G, **kw)
            rob.append({"scenario": label, "G_pct": G, "route": route, "yield_in_range": round(y, 4),
                        "breakeven_final": round(breakeven_y(yref, y, p, "final"), 3) if route != REF else 1.0,
                        "breakeven_probe": round(breakeven_y(yref, y, p, "probe"), 3) if route != REF else 1.0})
with open(os.path.join(DATA, "process_robustness.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rob[0].keys())); w.writeheader(); w.writerows(rob)

json.dump(summary, open(os.path.join(DATA, "cost_summary.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
