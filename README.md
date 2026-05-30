# Rocket Propellant Combustion Analysis Tool

## Live demo

🚀 **[Open interactive dashboard](https://mkws2026-rocket-propellant.streamlit.app)**

Course project for **Metody komputerowe w spalaniu** (Computer Methods in
Combustion) — Warsaw University of Technology, Faculty of Power and
Aeronautical Engineering.

A Cantera-based engine that computes chamber conditions (constant-enthalpy
equilibrium) and nozzle exit conditions (isentropic expansion with frozen
entropy + equilibrium at the exit) for five rocket propellant combinations.
Results are exported to CSV, plotted as six publication-quality matplotlib
figures, and made explorable through an interactive Streamlit dashboard. The
calculated specific impulses are validated against five real engines
(SpaceX Raptor, RS-25, Merlin 1D, RD-180, Vulcain 2).

## Installation

Requires **Python 3.10+**.

```bash
pip install -r requirements.txt
```

## How to run

**1. Generate the CSV results and all six figures**

```bash
python propellant_analysis.py     # writes results/*.csv  (~1–2 min)
python plots.py                   # writes figures/*.png
```

`python plots.py` will run the analysis automatically if the CSVs are missing,
so a fresh checkout only needs the second command.

**2. Launch the interactive dashboard**

```bash
streamlit run app.py
```

Then open the URL Streamlit prints (usually <http://localhost:8501>).

## Files

| File | What it does |
|---|---|
| `propellant_analysis.py` | Core Cantera engine: chamber HP equilibrium, isentropic SP nozzle expansion, parameter sweeps, CSV export. |
| `plots.py` | Six matplotlib figures at 300 dpi: Isp vs O/F, T vs O/F, species composition bars, (O/F, p) heatmap, density-Isp chart, c* vs O/F. |
| `app.py` | Streamlit dashboard with sliders for propellant, O/F, chamber pressure, and exit pressure; live metrics, Plotly charts, and a real-engine validation table. |
| `requirements.txt` | Python dependencies. |
| `results/` | Auto-generated CSVs (`results_p100bar.csv`, `heatmap_data.csv`). |
| `figures/` | Auto-generated PNG figures. |

## Propellants analysed

All chemistry uses **GRI-Mech 3.0** (`gri30.yaml`). Approximate peak values are
at p_c = 100 bar with isentropic expansion to vacuum.

| Propellant     | Fuel  | Oxidiser | O/F range | Peak Isp_vac [s] | Peak T_chamber [K] |
|----------------|-------|----------|-----------|------------------|---------------------|
| LOX/LH2        | H2    | O2       | 2.0–10.0  | ~510             | ~3650               |
| LOX/CH4        | CH4   | O2       | 1.5–6.0   | ~430             | ~3690               |
| LOX/C2H6       | C2H6  | O2       | 1.0–4.0   | ~425             | ~3750               |
| N2O4/NH3       | NH3   | NO2      | 0.5–3.0   | ~370             | ~3200               |
| N2O/C2H4       | C2H4  | N2O      | 2.0–8.0   | ~350             | ~3560               |

### Species substitutions

GRI-Mech 3.0 was developed for natural-gas combustion and does not contain
ethanol, hydrazine, or kerosene. The following substitutions were made:

- **C2H5OH (ethanol) → C2H6 (ethane)** — same carbon count, closest hydrocarbon
  analogue. Used in place of the spec's LOX/C2H5OH.
- **N2H4 (hydrazine) → NH3 (ammonia)** — N2H4 decomposes to NH3 + N2 + H2 in
  catalytic thrusters; NH3 captures the dominant chemistry.
- **N2O4 (dinitrogen tetroxide) → NO2** — N2O4 ⇌ 2 NO2 dissociates almost
  completely at chamber temperatures, so NO2 is the correct species to put in.
- **RP-1 (kerosene) → C2H6** — used only for the real-engine validation table
  (Merlin 1D, RD-180). Real RP-1 has a higher boiling point and slightly higher
  density Isp than ethane, so the calculated values are an approximation.

### Why the calculated Isp is higher than real engines

This is an **idealised equilibrium model**. Real engines lose Isp to:

- finite chemical kinetics (recombination during nozzle expansion is incomplete)
- boundary-layer and film-cooling losses
- finite nozzle expansion (real area ratios are finite, not infinite)
- combustion inefficiency in the injector / chamber

A 5–20 % overshoot vs. real-engine Isp is therefore expected and is visible in
the validation table inside the Streamlit app.
