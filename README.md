# Cost-effective fabrication of MEMS pressure sensors

Code and data for my research paper on how the choice of membrane fabrication route changes the cost of
one calibrated MEMS pressure sensor (Lumiere Education research program, 2026).

The paper follows one chain: membrane thickness → sensitivity → share of chips the readout chip can
still correct → cost per calibrated sensor. The reference device is the silicon-on-nothing sensor of
Su et al. (2019): 600 × 600 µm membrane, 13 µm thick, 0–100 kPa.

## What is in the repository

| Path | What it is |
| --- | --- |
| `membrane_fem.py` | 3D finite-element model of the membrane in scikit-fem (anisotropic silicon, quadratic hexahedra). One run = one command. |
| `nonlinear_check.py` | large-deflection version of the clamped model, used only as a check |
| `process_cost.py` | Monte Carlo of thickness spread for each route and the cost model; reads `data/public/` |
| `collect.py` | turns the raw run files in `runs/` into the CSV tables in `data/` |
| `make_figures.py` | figures 1–7 |
| `process_flow.py` | cross-sections of the SON fabrication flow (figure 10) |
| `fab_figures.py` | chip layout and back-end flow (figures 11 and 12) |
| `fem3d/` | independent check in CalculiX with a Gmsh mesh (tetrahedra), 3D viewer and figures 8–9 |
| `runs/` | argument lists and raw output of every scikit-fem run |
| `data/` | results as CSV files |
| `data/public/` | input data from public sources, with a table of how each source fits this study |
| `results/fem3d/` | CalculiX results, a ParaView file (`.vtu`) and an interactive 3D page (`.html`) |
| `figures/` | all figures |
| `notebook/` | Google Colab notebook that re-runs the main numbers in a browser |

## Running it

Python 3.10 or newer.

```
pip install -r requirements.txt
python membrane_fem.py --iso 169 0.3 --nxy 12 --nz 2
python membrane_fem.py --nxy 12 --nz 2 --frame 150 100
python process_cost.py
python collect.py
python make_figures.py
python process_flow.py
python fab_figures.py
```

The first command is the solver check against plate theory (about 0.485 µm). The second is the
13 µm reference membrane on a silicon frame (0.6277 µm, 34.56 mV/V at 100 kPa). All lengths are in µm and
pressures in kPa; `python membrane_fem.py --help` lists the options. Each full list of runs can be repeated
with `runs/queue.sh listA.txt resA.txt` and so on.

The CalculiX check needs `ccx` (version 2.17 was used) on the PATH:

```
cd fem3d
./run_ccx.sh
python compare.py
python view3d.py runs/fine6coarse30frame150100
```

The Monte Carlo uses fixed random seeds, so the numbers come out the same on every computer. I ran the
scikit-fem cases on a cloud machine and on my laptop and got identical results, and the CalculiX runs on
my laptop.

## Main numbers

| Check | Result |
| --- | --- |
| isotropic plate vs plate theory | 0.485 vs 0.480 µm |
| reference membrane vs ANSYS (Su et al. 2019) | 0.628 vs 0.64 µm |
| scikit-fem vs CalculiX, 9–17 µm | output within 0.6 %, deflection within 0.5 % |
| sensitivity vs thickness | power law, exponent −2.07 |
| timed KOH etch, chips within ±10 % | 79 % |
| break-even wafer price, ±10 %, final-test screening | 2.6 × timed-etch wafer |

## License

Code: MIT (see `LICENSE`). Data and figures: CC BY 4.0. If you use them, please cite the repository
(`CITATION.cff`).
