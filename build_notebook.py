import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
files = ["membrane_fem.py", "nonlinear_check.py", "process_cost.py", "data/sweep_su.csv"]
files += ["data/public/" + f for f in sorted(os.listdir(os.path.join(HERE, "data", "public"))) if f.endswith(".csv")]
blob = {f: open(os.path.join(HERE, f), encoding="utf-8").read() for f in files}

setup = "import os\n!pip -q install scikit-fem==12.0.2\nFILES = " + repr(blob) + "\n"
setup += ("for name, text in FILES.items():\n"
          "    os.makedirs(os.path.dirname(name) or '.', exist_ok=True)\n"
          "    open(name, 'w', encoding='utf-8').write(text)\n"
          "print(len(FILES), 'files written')")

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
nb = nbf.v4.new_notebook()
nb.cells = [
    md("# MEMS pressure sensor: re-run the main results\n"
       "Upload to Google Colab and use Runtime > Run all. The first cell writes the model scripts and the input tables "
       "(same files as in the repository), the next cells re-run the checks. Expected numbers are written above each cell."),
    code(setup),
    md("**Solver check**, isotropic plate E = 169 GPa, ν = 0.3, 600 µm, 13 µm, 100 kPa. Plate theory gives 0.480 µm, "
       "the paper reports w_centre_um = 0.4847."),
    code("!python membrane_fem.py --iso 169 0.3 --nxy 12 --nz 2"),
    md("**Reference membrane**, anisotropic Si, 13 µm, 100 kPa. Expected: clamp 0.5743 µm and 36.34 mV/V, "
       "frame 0.6277 µm and 34.56 mV/V (Su et al. 2019 report 0.64 µm from ANSYS)."),
    code("!python membrane_fem.py --nxy 12 --nz 2\n!python membrane_fem.py --nxy 12 --nz 2 --frame 150 100"),
    md("**Singh et al. (2014) 12 µm membrane**, 1097 µm wide, 50 kPa, stress averaged over the 1 µm resistor depth. "
       "Expected 34.89 mV/V, which is 69.8 mV/V/bar. Their FEM 70.5, measured 71.9–84.5 mV/V/bar. Takes 2–4 min."),
    code("!python membrane_fem.py --nxy 12 --nz 2 --t 12 --side 1097 --p 50 --pi 3 -1 80 "
         "--resL 60 10 10 --resT 60 10 20 --frame 300 600 --depth 1"),
    md("**Large deflection**, same coarse mesh for both. Expected linear 0.5702 µm / 36.33 mV/V, "
       "nonlinear 0.5697 µm / 36.52 mV/V."),
    code("!python membrane_fem.py --nxy 8 --nz 1\n!python nonlinear_check.py --nxy 8 --nz 1 --steps 1"),
    md("**Monte Carlo and cost model**. Expected at ±10 %: timed etch keeps 79 % of chips, break-even 2.6× "
       "(final test) and 1.35× (wafer probe)."),
    code("!python process_cost.py\nimport pandas as pd\n"
         "display(pd.read_csv('data/process_montecarlo.csv'))\n"
         "display(pd.read_csv('data/cost_breakeven.csv').query('G_pct == 10'))\n"
         "display(pd.read_csv('data/process_robustness.csv').query('G_pct == 10'))"),
]
os.makedirs(os.path.join(HERE, "notebook"), exist_ok=True)
nbf.write(nb, os.path.join(HERE, "notebook", "MEMS_membrane_check.ipynb"))
print("notebook/MEMS_membrane_check.ipynb")
