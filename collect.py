import json, glob, csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data"); os.makedirs(DATA, exist_ok=True)

runs = []; nl_runs = []
for f in sorted(glob.glob(os.path.join(HERE, "res*.txt")) + glob.glob(os.path.join(HERE, "runs", "res*.txt"))):
    for line in open(f):
        if "|" not in line: continue
        a, b = line.split("|", 1)
        try: d = json.loads(b)
        except Exception: continue
        d["args"] = a.strip(); d["batch"] = os.path.basename(f)
        if d.get("model", "").startswith("nonlinear"): nl_runs.append(d); continue
        runs.append(d)

def bc(d): return "frame" if d.get("frame_um") else "clamp"
cols = ["batch", "args", "material", "side_um", "t_um", "p_kPa", "nxy", "nz", "dofs", "frame_um", "w_centre_um", "w_over_t",
        "sigma_edge_MPa", "sigma_x_2um_in_MPa", "sigma_centre_MPa", "sxL_MPa", "syL_MPa", "sxT_MPa", "syT_MPa",
        "res_depth_um", "dR_L", "dR_T", "vout_mV_per_V", "sens_mV_per_V_per_kPa", "pi_1e11", "resL_um", "resT_um", "seconds"]
with open(os.path.join(DATA, "all_fem_runs.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
    for d in runs: w.writerow({k: d.get(k) for k in cols})

def is_su(d): return d["side_um"] == 600.0 and d["p_kPa"] == 100.0 and d["pi_1e11"] == [6.6, -1.1, 138.1] and d["material"].startswith("aniso")
for name, want in (("sweep_su.csv", "frame"), ("sweep_su_clamp.csv", "clamp")):
    rows = {}
    for d in runs:
        if is_su(d) and d["nxy"] == 12 and d["nz"] == 2 and bc(d) == want and not d.get("res_depth_um"):
            if want == "frame" and d["frame_um"] != [150.0, 100.0]: continue
            rows[d["t_um"]] = d
    with open(os.path.join(DATA, name), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["t_um", "w_centre_um", "sigma_edge_MPa", "sxL_MPa", "vout_mV_per_V", "sens_mV_per_V_per_kPa"])
        for t in sorted(rows):
            d = rows[t]; w.writerow([t, round(d["w_centre_um"], 5), round(d["sigma_edge_MPa"], 3), round(d["sxL_MPa"], 3),
                                     round(d["vout_mV_per_V"], 4), round(d["sens_mV_per_V_per_kPa"], 6)])
    print(name, len(rows), "points")

singh_sim = {9: 106.68, 10: 96.72, 11: 82.52, 12: 70.49, 13: 60.07, 14: 51.68, 15: 44.82}
singh_meas = {10: (86.7, 96.6), 12: (71.9, 84.5), 14: (50.5, 54.5)}
with open(os.path.join(DATA, "singh_comparison.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["t_um", "side_um", "edge_condition", "nxy", "res_depth_um", "this_model_mV_per_V_per_bar", "singh_simulated", "singh_measured_low", "singh_measured_high", "w_over_t"])
    for d in sorted(runs, key=lambda d: (d["t_um"], bc(d), d["nxy"])):
        if d["pi_1e11"] == [3.0, -1.0, 80.0]:
            t = int(round(d["t_um"])); s = d["vout_mV_per_V"] / (d["p_kPa"] / 100.0)
            m = singh_meas.get(t, ("", ""))
            w.writerow([d["t_um"], d["side_um"], bc(d), d["nxy"], d.get("res_depth_um", 0.0) or 0.0, round(s, 2), singh_sim.get(t, ""), m[0], m[1], round(d["w_over_t"], 3)])
print("runs", len(runs))

lin = [d for d in runs if d["batch"] == "resNL.txt"]
with open(os.path.join(DATA, "nonlinear_check.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["case", "t_um", "side_um", "p_kPa", "w_linear_um", "w_nonlinear_um", "w_change_pct",
                                    "vout_linear_mV_per_V", "vout_nonlinear_mV_per_V", "vout_change_pct", "w_over_t_linear"])
    for n in nl_runs:
        l = next(d for d in lin if d["t_um"] == n["t_um"] and d["side_um"] == n["side_um"])
        case = "Su membrane" if n["side_um"] == 600.0 else "Singh membrane"
        w.writerow([case, n["t_um"], n["side_um"], n["p_kPa"], round(l["w_centre_um"], 4), round(n["w_centre_um"], 4),
                    round(100 * (n["w_centre_um"] / l["w_centre_um"] - 1), 1), round(l["vout_mV_per_V"], 3),
                    round(n["vout_mV_per_V"], 3), round(100 * (n["vout_mV_per_V"] / l["vout_mV_per_V"] - 1), 1), round(l["w_over_t"], 3)])
print("nonlinear pairs", len(nl_runs))
