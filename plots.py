"""
Publication-quality plots for the rocket propellant analysis.

Reads CSVs produced by propellant_analysis.py and writes six PNG figures
at 300 dpi into the figures/ directory. If the CSVs are missing it runs
the analysis first.
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import propellant_analysis as pa


FIGURES_DIR = 'figures'
MAIN_CSV    = os.path.join(pa.RESULTS_DIR, 'results_p100bar.csv')
HEATMAP_CSV = os.path.join(pa.RESULTS_DIR, 'heatmap_data.csv')


COLORS = {
    'LOX/LH2':  'tab:blue',
    'LOX/CH4':  'tab:orange',
    'LOX/C2H6': 'tab:green',
    'N2O4/NH3': 'tab:red',
    'N2O/C2H4': 'tab:purple',
}


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

def setup_style() -> None:
    for style in ('seaborn-v0_8-paper', 'seaborn-paper', 'seaborn-v0_8'):
        try:
            plt.style.use(style)
            break
        except OSError:
            continue
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'legend.fontsize': 9,
        'mathtext.fontset': 'cm',
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'figure.dpi': 110,
    })


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not (os.path.exists(MAIN_CSV) and os.path.exists(HEATMAP_CSV)):
        print('CSVs missing — running propellant_analysis.main() first.')
        pa.main()
    return pd.read_csv(MAIN_CSV), pd.read_csv(HEATMAP_CSV)


def _color(name: str) -> str:
    return COLORS.get(name, 'gray')


# ---------------------------------------------------------------------------
# Plot 1 — Isp_vac vs O/F
# ---------------------------------------------------------------------------

def plot_isp_vs_of(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, sub in df.groupby('propellant_name', sort=False):
        sub = sub.sort_values('OF_ratio')
        ax.plot(sub['OF_ratio'], sub['Isp_vac'], label=name,
                color=_color(name), linewidth=2)
        # stoichiometric line for this propellant
        p = next(p for p in pa.PROPELLANTS if p['name'] == name)
        of_stoich = pa.stoichiometric_OF(p['fuel'], p['oxidizer'])
        if sub['OF_ratio'].min() <= of_stoich <= sub['OF_ratio'].max():
            ax.axvline(of_stoich, color=_color(name), linestyle=':',
                       linewidth=0.8, alpha=0.6)
    # real engines as star markers
    for eng_name, eng in pa.REAL_ENGINES.items():
        ax.plot(eng['OF'], eng['Isp_vac'], marker='*', markersize=14,
                color='black', markeredgecolor='white', markeredgewidth=0.7,
                linestyle='none')
        ax.annotate(eng_name, (eng['OF'], eng['Isp_vac']),
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    ax.set_xlabel('O/F ratio [-]')
    ax.set_ylabel(r'$I_{sp,vac}$  [s]')
    ax.set_title('Vacuum specific impulse vs O/F ratio  ($p_c = 100$ bar)')
    ax.legend(loc='lower right', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    out = os.path.join(FIGURES_DIR, 'fig1_isp_vs_of.png')
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f'wrote {out}')


# ---------------------------------------------------------------------------
# Plot 2 — Adiabatic flame temperature vs O/F
# ---------------------------------------------------------------------------

def plot_temperature_vs_of(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, sub in df.groupby('propellant_name', sort=False):
        sub = sub.sort_values('OF_ratio')
        ax.plot(sub['OF_ratio'], sub['T_chamber'], label=name,
                color=_color(name), linewidth=2)
    ax.axhline(3000, color='gray', linestyle='--', linewidth=0.8)
    ax.text(ax.get_xlim()[1] * 0.98, 3010, 'typical material limit 3000 K',
            fontsize=8, color='gray', ha='right')
    ax.axhline(3700, color='red', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.text(ax.get_xlim()[1] * 0.98, 3710, 'regen-cooled limit ~3700 K',
            fontsize=8, color='red', ha='right', alpha=0.7)
    ax.set_xlabel('O/F ratio [-]')
    ax.set_ylabel(r'$T_{chamber}$  [K]')
    ax.set_title('Adiabatic flame temperature vs O/F ratio  ($p_c = 100$ bar)')
    ax.legend(loc='lower right', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    out = os.path.join(FIGURES_DIR, 'fig2_temperature_vs_of.png')
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f'wrote {out}')


# ---------------------------------------------------------------------------
# Plot 3 — Top 8 species at the optimal O/F (per propellant)
# ---------------------------------------------------------------------------

def plot_species_composition(df: pd.DataFrame) -> None:
    propellants = list(df['propellant_name'].unique())
    n = len(propellants)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 4.0 * nrows))
    axes = np.atleast_2d(axes).flatten()

    for ax, name in zip(axes, propellants):
        sub = df[df['propellant_name'] == name]
        i = sub['Isp_vac'].idxmax()
        row = sub.loc[i]
        names = [row[f'species_{k}_name'] for k in range(1, pa.N_TOP_SPECIES + 1)]
        fracs = [row[f'species_{k}_frac'] for k in range(1, pa.N_TOP_SPECIES + 1)]
        pairs = [(n2, f) for n2, f in zip(names, fracs) if n2 and f > 0]
        pairs.sort(key=lambda kv: kv[1], reverse=True)
        sn = [k for k, _ in pairs]
        sf = [v for _, v in pairs]
        bars = ax.bar(range(len(sn)), sf, color=_color(name), alpha=0.85,
                      edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(sn)))
        ax.set_xticklabels(sn, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel('mole fraction [-]')
        ax.set_title(f"{name}\nOF*={row['OF_ratio']:.2f},  T={row['T_chamber']:.0f} K",
                     fontsize=10)
        ax.set_ylim(0, max(sf) * 1.18)
        ax.grid(True, axis='y', alpha=0.3)
        for b, v in zip(bars, sf):
            ax.text(b.get_x() + b.get_width() / 2, v, f'{v:.2f}',
                    ha='center', va='bottom', fontsize=7)
    for ax in axes[len(propellants):]:
        ax.axis('off')
    fig.suptitle('Chamber composition (top 8 species) at optimal O/F',
                 fontsize=13, y=1.00)
    fig.tight_layout()
    out = os.path.join(FIGURES_DIR, 'fig3_species_composition.png')
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f'wrote {out}')


# ---------------------------------------------------------------------------
# Plot 4 — Heatmap of Isp_vac (O/F vs p_chamber) for LOX/LH2
# ---------------------------------------------------------------------------

def plot_heatmap_isp(heatmap_df: pd.DataFrame, propellant_name: str = 'LOX/LH2') -> None:
    sub = heatmap_df[heatmap_df['propellant_name'] == propellant_name]
    if sub.empty:
        print(f'heatmap: no data for {propellant_name}; skipping')
        return
    pivot = sub.pivot_table(index='p_chamber_bar', columns='OF_ratio',
                            values='Isp_vac')
    pivot = pivot.sort_index().sort_index(axis=1)
    of_vals = pivot.columns.values.astype(float)
    p_vals  = pivot.index.values.astype(float)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    mesh = ax.pcolormesh(of_vals, p_vals, pivot.values,
                         shading='auto', cmap='RdYlGn')
    # contour lines at notable Isp values
    levels = [400, 420, 440, 460, 480]
    levels = [L for L in levels
              if np.nanmin(pivot.values) <= L <= np.nanmax(pivot.values)]
    if levels:
        cs = ax.contour(of_vals, p_vals, pivot.values, levels=levels,
                        colors='black', linewidths=0.8)
        ax.clabel(cs, inline=True, fmt='%d s', fontsize=8)

    # real engines that use this propellant
    for eng_name, eng in pa.REAL_ENGINES.items():
        if eng['sim_propellant'] == propellant_name:
            ax.plot(eng['OF'], eng['p_bar'], marker='*', markersize=18,
                    color='white', markeredgecolor='black', markeredgewidth=1.0,
                    linestyle='none')
            ax.annotate(eng_name, (eng['OF'], eng['p_bar']),
                        xytext=(8, 4), textcoords='offset points',
                        fontsize=9, color='black',
                        bbox=dict(boxstyle='round,pad=0.2',
                                  facecolor='white', alpha=0.8, edgecolor='none'))

    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(r'$I_{sp,vac}$  [s]')
    ax.set_xlabel('O/F ratio [-]')
    ax.set_ylabel(r'$p_{chamber}$  [bar]')
    ax.set_title(f'{propellant_name}: vacuum $I_{{sp}}$ vs (O/F, $p_c$)')
    out = os.path.join(FIGURES_DIR, 'fig4_heatmap_isp.png')
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f'wrote {out}')


# ---------------------------------------------------------------------------
# Plot 5 — Density-Isp vs Isp (propellant selection chart)
# ---------------------------------------------------------------------------

def plot_density_isp_chart(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for name, sub in df.groupby('propellant_name', sort=False):
        i = sub['Isp_vac'].idxmax()
        row = sub.loc[i]
        # rho_isp_vac is kg*s/m^3 in CSV; divide by 1000 -> kg*s/L
        rho_isp_kgL = row['rho_isp_vac'] / 1000.0
        ax.scatter(row['Isp_vac'], rho_isp_kgL,
                   color=_color(name), s=160, edgecolor='black', linewidth=0.7,
                   zorder=3, label=name)
        ax.annotate(name, (row['Isp_vac'], rho_isp_kgL),
                    xytext=(8, 6), textcoords='offset points', fontsize=10)

    ax.set_xlabel(r'$I_{sp,vac}$  [s]    (efficiency — how far per kg of propellant)')
    ax.set_ylabel(r'$\rho\,I_{sp}$  [kg$\cdot$s / L]    (compactness — thrust per tank volume)')
    ax.set_title('Density-$I_{sp}$ vs $I_{sp,vac}$  —  propellant selection chart')
    ax.grid(True, alpha=0.3)
    ax.text(0.02, 0.98,
            'upper-right  =  high performance AND volumetrically dense\n'
            'lower-right  =  high $I_{sp}$ but bulky tanks (LH$_2$ trade-off)',
            transform=ax.transAxes, fontsize=8.5, va='top',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                      edgecolor='gray', alpha=0.85))
    out = os.path.join(FIGURES_DIR, 'fig5_density_isp_chart.png')
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f'wrote {out}')


# ---------------------------------------------------------------------------
# Plot 6 — c* (characteristic velocity) vs O/F
# ---------------------------------------------------------------------------

def plot_cstar_vs_of(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, sub in df.groupby('propellant_name', sort=False):
        sub = sub.sort_values('OF_ratio')
        ax.plot(sub['OF_ratio'], sub['c_star'], label=name,
                color=_color(name), linewidth=2)
    ax.set_xlabel('O/F ratio [-]')
    ax.set_ylabel(r'$c^*$  [m/s]')
    ax.set_title('Characteristic velocity $c^*$ vs O/F ratio  ($p_c = 100$ bar)')
    ax.legend(loc='lower right', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    out = os.path.join(FIGURES_DIR, 'fig6_cstar_vs_of.png')
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f'wrote {out}')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    os.makedirs(FIGURES_DIR, exist_ok=True)
    setup_style()
    df, heatmap = load_data()
    plot_isp_vs_of(df)
    plot_temperature_vs_of(df)
    plot_species_composition(df)
    plot_heatmap_isp(heatmap, 'LOX/LH2')
    plot_density_isp_chart(df)
    plot_cstar_vs_of(df)
    print(f'\nAll figures saved to {FIGURES_DIR}/')


if __name__ == '__main__':
    main()
