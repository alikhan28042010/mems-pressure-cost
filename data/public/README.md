# Input data from public sources

No foundry gave me price quotes, so the cost model uses numbers that anyone can look up. Every file here
has a source column. `process_cost.py` reads `cost_model_inputs.csv`, `thickness_control.csv` and
`thickness_scenarios.csv` directly; the other files are used for checks in the paper.

| File | What it holds | Used for |
| --- | --- | --- |
| `cost_model_inputs.csv` | wafer, package, readout chip, calibration and yield inputs with low/base/high values | cost model |
| `thickness_control.csv` | membrane thickness spread of each fabrication route, the die size each route needs | Monte Carlo and cost model |
| `thickness_scenarios.csv` | spreads tried for APSM and MIS, which have no published data | Monte Carlo |
| `validation_points.csv` | measured and simulated values from Su et al. (2019), Singh et al. (2014) and Kloeck et al. (1989) | model checks |
| `lawes_2007_values.csv` | numbers from the Lawes (2007) MEMS cost model | where the wafer price comes from |
| `mpw_prices_2026.csv` | EUROPRACTICE multi-project wafer prices, 2026 | prototype cost, for comparison only |
| `market_prices.csv` | Mouser prices of the Bosch BMP390 barometric sensor, 31 Oct 2026 | upper limit for the cost of a finished sensor |

## Do the sources fit this study?

The study is one sensor: a piezoresistive absolute pressure sensor for air, 0–100 kPa, 0–80 °C, with a
600 × 600 µm silicon membrane 13 µm thick, 1.2 mm die (2.0 mm for back-side KOH routes) on a 150 mm wafer.

| Source | Their case | This study | Same? |
| --- | --- | --- | --- |
| Lawes (2007), $800 per wafer | CMOS + MEMS accelerometer, 150 mm wafer, 10,000 wafers per year | piezoresistive pressure sensor, 150 mm wafer | wafer size and volume match; device and year do not, so the range is $400–1,600 |
| Fitzgerald (2021), package < $0.10 | molded package for a pressure sensor | pressure sensor | yes |
| Pavlin & Novak (2008), 72 min per 128 sensors | piezoresistive pressure sensors, 3 pressures × 3 temperatures | piezoresistive, 0–80 °C needs several temperatures | yes |
| Singh et al. (2014), ±1 µm | KOH-etched membranes 10–14 µm thick | 9–17 µm, nominal 13 µm | yes |
| Kloeck et al. (1989), 0.2 µm s.d. | electrochemical etch stop on epitaxial silicon membranes | same route | yes |
| Okmetic, ±0.5 / ±0.2 µm | cavity-SOI device layer | same route | yes |
| Su et al. (2019), ±1 µm and < ±0.2 µm | cavity-SOI and SON membranes, 13 µm | the reference sensor itself | yes |
| EUROPRACTICE 2026 MPW prices | prototype runs, a few dies per design | volume production | no, only used to show the gap between prototype and volume cost |
| Mouser BMP390 price | finished piezoresistive barometer with ASIC, 2 × 2 mm package | same application | it includes the maker's margin, so it is only an upper limit |

What is still missing: the price of cavity-SOI wafers (the paper reports a break-even price instead),
thickness data for APSM and MIS, and real foundry quotes.

The data in this folder are shared under CC BY 4.0.
