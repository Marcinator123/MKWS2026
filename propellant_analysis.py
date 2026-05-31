"""
Rocket Propellant Combustion Analysis — Core Cantera Engine
-----------------------------------------------------------
Computes chamber conditions (HP equilibrium) and nozzle exit conditions
(isentropic + SP equilibrium) for several rocket propellant combinations
using the GRI-Mech 3.0 mechanism.

Species substitutions (GRI-Mech 3.0 does not contain ethanol or hydrazine):
    C2H5OH (ethanol)  ->  C2H6  (ethane)   — closest hydrocarbon analogue
    N2H4   (hydrazine) -> NH3   (ammonia)  — main decomposition product
    RP-1   (kerosene) ->  C2H6  (ethane)   — for real-engine validation only

Course: Metody komputerowe w spalaniu (Computer Methods in Combustion)
Warsaw University of Technology, Faculty of Power and Aeronautical Engineering
"""

from __future__ import annotations

import math
import os
import warnings
from typing import Iterable

import cantera as ct
import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **_kwargs):
        return x


# ---------------------------------------------------------------------------
# Physical constants & configuration
# ---------------------------------------------------------------------------

G0 = 9.80665                # standard gravity [m/s^2]
T0 = 300.0                  # initial reactants temperature [K]
MECHANISM = 'gri30.yaml'
RESULTS_DIR = 'results'

# How many species columns to store per row (used by plots.py for composition bars)
N_TOP_SPECIES = 8


# ---------------------------------------------------------------------------
# Propellant registry
# ---------------------------------------------------------------------------

PROPELLANTS = [
    {'name': 'LOX/LH2',   'fuel': 'H2',   'oxidizer': 'O2',
     'OF_range': (2.0, 10.0), 'note': ''},
    {'name': 'LOX/CH4',   'fuel': 'CH4',  'oxidizer': 'O2',
     'OF_range': (1.5, 6.0),  'note': ''},
    {'name': 'LOX/C2H6',  'fuel': 'C2H6', 'oxidizer': 'O2',
     'OF_range': (1.0, 4.0),
     'note': 'ethane substituted for C2H5OH (ethanol not in gri30)'},
    {'name': 'N2O4/NH3',  'fuel': 'NH3',  'oxidizer': 'NO2',
     'OF_range': (0.5, 3.0),
     'note': 'NH3 substituted for N2H4; NO2 used as N2O4 dissociation product'},
    {'name': 'N2O/C2H4',  'fuel': 'C2H4', 'oxidizer': 'N2O',
     'OF_range': (0.5, 10.0),  'note': ''},
]


# Real engine reference data for validation.
# 'sim_propellant' is the propellant in our PROPELLANTS list used to simulate
# this engine (RP-1 is mapped to LOX/C2H6 because RP-1 is not in gri30).
REAL_ENGINES = {
    'SpaceX Raptor':   {'propellant': 'LOX/CH4',  'sim_propellant': 'LOX/CH4',
                        'OF': 3.55, 'Isp_vac': 363, 'p_bar': 300},
    'RS-25 (Shuttle)': {'propellant': 'LOX/LH2',  'sim_propellant': 'LOX/LH2',
                        'OF': 6.03, 'Isp_vac': 453, 'p_bar': 207},
    'Merlin 1D':       {'propellant': 'LOX/RP-1', 'sim_propellant': 'LOX/C2H6',
                        'OF': 2.36, 'Isp_vac': 348, 'p_bar': 97},
    'RD-180':          {'propellant': 'LOX/RP-1', 'sim_propellant': 'LOX/C2H6',
                        'OF': 2.72, 'Isp_vac': 338, 'p_bar': 267},
    'Vulcain 2':       {'propellant': 'LOX/LH2',  'sim_propellant': 'LOX/LH2',
                        'OF': 6.10, 'Isp_vac': 434, 'p_bar': 116},
}


# Exit pressures used in the main sweep
P_VAC_BAR = 0.001    # vacuum (space)
P_SL_BAR  = 1.013    # sea level (Earth surface)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _new_gas() -> ct.Solution:
    return ct.Solution(MECHANISM)


def species_available(gas: ct.Solution, name: str) -> bool:
    return name in gas.species_names


def filter_available_propellants(propellants: Iterable[dict]) -> list[dict]:
    """Return propellants whose fuel and oxidizer are both in the mechanism."""
    gas = _new_gas()
    kept = []
    for p in propellants:
        missing = [s for s in (p['fuel'], p['oxidizer']) if not species_available(gas, s)]
        if missing:
            print(f"[skip propellant] {p['name']}: species {missing} not in {MECHANISM}")
        else:
            kept.append(p)
    return kept


def stoichiometric_OF(fuel: str, oxidizer: str) -> float:
    """Analytical stoichiometric oxidizer-to-fuel mass ratio via Cantera."""
    gas = _new_gas()
    gas.TP = T0, ct.one_atm
    gas.set_equivalence_ratio(1.0, fuel, oxidizer)
    y = gas.mass_fraction_dict()
    y_fuel = y.get(fuel, 0.0)
    y_ox   = y.get(oxidizer, 0.0)
    if y_fuel <= 0.0:
        return float('nan')
    return y_ox / y_fuel


def _top_n_species(gas: ct.Solution, n: int) -> list[tuple[str, float]]:
    items = sorted(gas.mole_fraction_dict().items(), key=lambda kv: kv[1], reverse=True)
    return items[:n]


# ---------------------------------------------------------------------------
# Step 1: chamber (HP equilibrium)
# ---------------------------------------------------------------------------

def compute_chamber(gas: ct.Solution, fuel: str, oxidizer: str,
                    OF: float, p_chamber_pa: float) -> dict:
    """Equilibrate at constant enthalpy & pressure; return chamber properties."""
    # mass fractions from O/F ratio
    y_fuel = 1.0 / (1.0 + OF)
    y_ox   = OF  / (1.0 + OF)
    gas.Y = {fuel: y_fuel, oxidizer: y_ox}
    gas.TP = T0, p_chamber_pa
    # ChemEquil prints noisy "Temperature outside valid range" warnings when
    # T exceeds 3000 K — equilibrium still converges, so suppress them.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        gas.equilibrate('HP')

    T_chamber   = gas.T
    gamma       = gas.cp_mass / gas.cv_mass
    M_mol       = gas.mean_molecular_weight                # kg/kmol
    R_specific  = gas.cp_mass - gas.cv_mass                # J/kg/K
    rho_chamber = gas.density_mass
    sound_speed = gas.sound_speed
    h0          = gas.enthalpy_mass
    s0          = gas.entropy_mass

    # Characteristic velocity c*  [m/s]
    Gamma_func = math.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    c_star = math.sqrt(gamma * R_specific * T_chamber) / Gamma_func

    top = _top_n_species(gas, N_TOP_SPECIES)

    return {
        'T_chamber': T_chamber,
        'gamma': gamma,
        'M_mol': M_mol,
        'R_specific': R_specific,
        'c_star': c_star,
        'rho_chamber': rho_chamber,
        'sound_speed': sound_speed,
        'h0': h0,
        's0': s0,
        'top_species': top,
    }


# ---------------------------------------------------------------------------
# Step 2: nozzle exit (isentropic + SP equilibrium)
# ---------------------------------------------------------------------------

def compute_exit(gas: ct.Solution, s0: float, h0: float,
                 p_exit_pa: float) -> dict:
    """Isentropic expansion to p_exit with chemical equilibrium at the exit."""
    try:
        gas.SP = s0, p_exit_pa
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            gas.equilibrate('SP')
        h_exit = gas.enthalpy_mass
        dh = h0 - h_exit
        if dh <= 0.0:
            warnings.warn(f'Non-physical expansion: h0={h0:.3e} <= h_exit={h_exit:.3e}')
            Ve = 0.0
        else:
            Ve = math.sqrt(2.0 * dh)
    except ct.CanteraError as e:
        warnings.warn(f'Exit equilibration failed at p_exit={p_exit_pa:.2f} Pa: {e}')
        Ve, h_exit = 0.0, float('nan')

    Isp = Ve / G0
    return {'Ve': Ve, 'Isp': Isp, 'h_exit': h_exit}


# ---------------------------------------------------------------------------
# Per-point analysis
# ---------------------------------------------------------------------------

def analyze_point(propellant: dict, OF: float, p_chamber_bar: float,
                  p_exit_vac_bar: float = P_VAC_BAR,
                  p_exit_sl_bar:  float = P_SL_BAR) -> dict:
    """Run a full chamber + exit calculation for one (propellant, OF, p) point."""
    gas = _new_gas()
    p_chamber_pa = p_chamber_bar * 1e5

    chamber = compute_chamber(gas, propellant['fuel'], propellant['oxidizer'],
                              OF, p_chamber_pa)

    # vacuum exit
    exit_vac = compute_exit(gas, chamber['s0'], chamber['h0'], p_exit_vac_bar * 1e5)
    # sea-level exit — fresh gas object because compute_exit mutates state
    gas_sl = _new_gas()
    chamber_sl = compute_chamber(gas_sl, propellant['fuel'], propellant['oxidizer'],
                                 OF, p_chamber_pa)
    exit_sl = compute_exit(gas_sl, chamber_sl['s0'], chamber_sl['h0'], p_exit_sl_bar * 1e5)

    # density-specific impulse [kg*s/m^3]; plots will convert to kg*s/L
    rho_isp_vac = chamber['rho_chamber'] * exit_vac['Isp']

    row = {
        'propellant_name': propellant['name'],
        'fuel': propellant['fuel'],
        'oxidizer': propellant['oxidizer'],
        'OF_ratio': OF,
        'p_chamber_bar': p_chamber_bar,
        'T_chamber': chamber['T_chamber'],
        'gamma': chamber['gamma'],
        'M_mol': chamber['M_mol'],
        'R_specific': chamber['R_specific'],
        'c_star': chamber['c_star'],
        'rho_chamber': chamber['rho_chamber'],
        'sound_speed': chamber['sound_speed'],
        'Ve_vac': exit_vac['Ve'],
        'Ve_sl':  exit_sl['Ve'],
        'Isp_vac': exit_vac['Isp'],
        'Isp_sl':  exit_sl['Isp'],
        'rho_isp_vac': rho_isp_vac,
    }
    for i, (name, frac) in enumerate(chamber['top_species'], start=1):
        row[f'species_{i}_name'] = name
        row[f'species_{i}_frac'] = frac
    # pad to N_TOP_SPECIES if equilibrium returned fewer significant species
    for i in range(len(chamber['top_species']) + 1, N_TOP_SPECIES + 1):
        row[f'species_{i}_name'] = ''
        row[f'species_{i}_frac'] = 0.0
    return row


# ---------------------------------------------------------------------------
# Sweeps
# ---------------------------------------------------------------------------

def run_main_analysis(p_chamber_bar: float = 100.0,
                      n_points: int = 40,
                      propellants: list[dict] | None = None) -> pd.DataFrame:
    """Sweep O/F across each propellant's range at fixed chamber pressure."""
    propellants = propellants or filter_available_propellants(PROPELLANTS)
    rows = []
    for p in propellants:
        of_min, of_max = p['OF_range']
        of_grid = np.linspace(of_min, of_max, n_points)
        for OF in tqdm(of_grid, desc=f"{p['name']:<10s} @ {p_chamber_bar:g} bar", leave=False):
            try:
                rows.append(analyze_point(p, float(OF), p_chamber_bar))
            except ct.CanteraError as e:
                print(f"[skip] {p['name']} OF={OF:.2f}: {e}")
    df = pd.DataFrame(rows)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out = os.path.join(RESULTS_DIR, f'results_p{int(p_chamber_bar)}bar.csv')
    df.to_csv(out, index=False)
    print(f'Wrote {out}  ({len(df)} rows)')
    return df


def run_heatmap_sweep(propellants: list[dict] | None = None,
                      pressures_bar: list[float] = (10, 30, 50, 100, 200, 300),
                      n_of: int = 25) -> pd.DataFrame:
    """Sweep both O/F and chamber pressure; minimal columns for heatmap plot."""
    propellants = propellants or filter_available_propellants(PROPELLANTS)
    rows = []
    for p in propellants:
        of_min, of_max = p['OF_range']
        of_grid = np.linspace(of_min, of_max, n_of)
        for p_bar in tqdm(pressures_bar, desc=f"heatmap {p['name']}", leave=False):
            for OF in of_grid:
                try:
                    r = analyze_point(p, float(OF), float(p_bar))
                    rows.append({
                        'propellant_name': p['name'],
                        'OF_ratio': r['OF_ratio'],
                        'p_chamber_bar': r['p_chamber_bar'],
                        'Isp_vac': r['Isp_vac'],
                        'T_chamber': r['T_chamber'],
                    })
                except ct.CanteraError as e:
                    print(f"[skip] {p['name']} OF={OF:.2f} p={p_bar}: {e}")
    df = pd.DataFrame(rows)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out = os.path.join(RESULTS_DIR, 'heatmap_data.csv')
    df.to_csv(out, index=False)
    print(f'Wrote {out}  ({len(df)} rows)')
    return df


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_summary(df: pd.DataFrame) -> None:
    print('\n' + '=' * 78)
    print('SUMMARY — peak vacuum Isp per propellant')
    print('=' * 78)
    if df.empty:
        print('(no data)')
        return
    rows = []
    for name, sub in df.groupby('propellant_name', sort=False):
        i = sub['Isp_vac'].idxmax()
        rows.append({
            'propellant': name,
            'OF*': round(sub.loc[i, 'OF_ratio'], 3),
            'Isp_vac [s]': round(sub.loc[i, 'Isp_vac'], 1),
            'Isp_sl [s]':  round(sub.loc[i, 'Isp_sl'], 1),
            'T_chamber [K]': round(sub.loc[i, 'T_chamber'], 0),
            'c* [m/s]': round(sub.loc[i, 'c_star'], 1),
        })
    print(pd.DataFrame(rows).to_string(index=False))
    print('=' * 78 + '\n')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> tuple[pd.DataFrame, pd.DataFrame]:
    print(f'Mechanism: {MECHANISM}')
    propellants = filter_available_propellants(PROPELLANTS)
    print(f'Analyzing {len(propellants)} propellants.\n')

    print('--- Main sweep @ p_chamber = 100 bar ---')
    df_main = run_main_analysis(p_chamber_bar=100.0, n_points=40,
                                propellants=propellants)
    print_summary(df_main)

    print('--- Heatmap sweep over (O/F, p_chamber) ---')
    df_heat = run_heatmap_sweep(propellants=propellants)

    print('\nStoichiometric O/F (analytical from gri30 molecular weights):')
    for p in propellants:
        print(f"  {p['name']:<10s}  OF_stoich = {stoichiometric_OF(p['fuel'], p['oxidizer']):.3f}")

    print('\nDone.')
    return df_main, df_heat


if __name__ == '__main__':
    main()
