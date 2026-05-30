"""
Rocket Propellant Combustion Analysis — Interactive Dashboard
=============================================================

To run:
    pip install -r requirements.txt
    streamlit run app.py

Course project — Metody komputerowe w spalaniu, Warsaw University of Technology.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import propellant_analysis as pa


# ---------------------------------------------------------------------------
# Page setup
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title='Rocket Propellant Analysis',
    layout='wide',
    page_icon=':rocket:',
)

st.title('Rocket Propellant Combustion Analysis')
st.caption('Cantera + GRI-Mech 3.0 — chamber equilibrium and isentropic '
           'nozzle expansion.  Course project, Warsaw University of Technology.')


PROPELLANT_NAMES = [p['name'] for p in pa.PROPELLANTS]
PROPELLANT_BY_NAME = {p['name']: p for p in pa.PROPELLANTS}


# ---------------------------------------------------------------------------
# Cached computations
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner='Computing point …')
def single_point(propellant_name: str, OF: float,
                 p_chamber_bar: float, p_exit_bar: float) -> dict:
    p = PROPELLANT_BY_NAME[propellant_name]
    # custom p_exit_vac = user value, p_exit_sl = 1.013 bar
    return pa.analyze_point(p, OF, p_chamber_bar,
                            p_exit_vac_bar=p_exit_bar,
                            p_exit_sl_bar=1.013)


@st.cache_data(show_spinner='Sweeping O/F …')
def of_sweep(propellant_name: str, p_chamber_bar: float,
             p_exit_bar: float, n: int = 30) -> pd.DataFrame:
    p = PROPELLANT_BY_NAME[propellant_name]
    of_min, of_max = p['OF_range']
    grid = np.linspace(of_min, of_max, n)
    rows = []
    for OF in grid:
        try:
            rows.append(pa.analyze_point(p, float(OF), p_chamber_bar,
                                         p_exit_vac_bar=p_exit_bar,
                                         p_exit_sl_bar=1.013))
        except Exception:
            continue
    return pd.DataFrame(rows)


@st.cache_data
def stoichiometric_point(propellant_name: str, p_chamber_bar: float,
                         p_exit_bar: float) -> dict:
    p = PROPELLANT_BY_NAME[propellant_name]
    OF = pa.stoichiometric_OF(p['fuel'], p['oxidizer'])
    return pa.analyze_point(p, OF, p_chamber_bar,
                            p_exit_vac_bar=p_exit_bar,
                            p_exit_sl_bar=1.013)


@st.cache_data(show_spinner='Validating against real engines …')
def real_engine_table(p_exit_bar: float = 0.001) -> pd.DataFrame:
    rows = []
    for eng_name, eng in pa.REAL_ENGINES.items():
        sim_name = eng['sim_propellant']
        if sim_name not in PROPELLANT_BY_NAME:
            continue
        try:
            r = pa.analyze_point(PROPELLANT_BY_NAME[sim_name],
                                 eng['OF'], eng['p_bar'],
                                 p_exit_vac_bar=p_exit_bar,
                                 p_exit_sl_bar=1.013)
            calc = r['Isp_vac']
            diff = (calc - eng['Isp_vac']) / eng['Isp_vac'] * 100.0
            rows.append({
                'Engine': eng_name,
                'Propellant (real)': eng['propellant'],
                'Propellant (sim)': sim_name,
                'O/F': eng['OF'],
                'p_c [bar]': eng['p_bar'],
                'Isp_vac real [s]': eng['Isp_vac'],
                'Isp_vac calc [s]': round(calc, 1),
                'Diff %': round(diff, 1),
            })
        except Exception as e:
            rows.append({'Engine': eng_name, 'error': str(e)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

st.sidebar.header('Inputs')

prop_name = st.sidebar.selectbox('Propellant', PROPELLANT_NAMES, index=0)
prop = PROPELLANT_BY_NAME[prop_name]
of_min, of_max = prop['OF_range']

of_default = float(pa.stoichiometric_OF(prop['fuel'], prop['oxidizer']))
of_default = float(np.clip(of_default, of_min, of_max))
OF = st.sidebar.slider('O/F ratio  [-]',
                       min_value=float(of_min), max_value=float(of_max),
                       value=of_default, step=0.1)

p_chamber = st.sidebar.slider('Chamber pressure  [bar]',
                              min_value=10, max_value=300, value=100, step=10)

p_exit = st.sidebar.slider('Exit pressure  [bar]',
                           min_value=0.001, max_value=2.0,
                           value=0.001, step=0.01)
st.sidebar.caption('0.001 bar = vacuum (space)   ·   1.013 bar = sea level')

if prop.get('note'):
    st.sidebar.info(f'**Substitution note**: {prop["note"]}')


# ---------------------------------------------------------------------------
# Compute current point + stoichiometric reference
# ---------------------------------------------------------------------------

current = single_point(prop_name, OF, p_chamber, p_exit)
stoich  = stoichiometric_point(prop_name, p_chamber, p_exit)


def delta_str(curr: float, ref: float, fmt: str = '{:+.1f}') -> str:
    return fmt.format(curr - ref)


# ---------------------------------------------------------------------------
# Metric cards
# ---------------------------------------------------------------------------

c1, c2, c3 = st.columns(3)

with c1:
    st.metric('Chamber temperature  [K]',
              f'{current["T_chamber"]:.0f}',
              delta_str(current['T_chamber'], stoich['T_chamber'], '{:+.0f}'))
    st.metric('Gamma  (cp/cv)  [-]',
              f'{current["gamma"]:.3f}',
              delta_str(current['gamma'], stoich['gamma'], '{:+.3f}'))

with c2:
    st.metric(f'Isp at p_exit = {p_exit:g} bar  [s]',
              f'{current["Isp_vac"]:.1f}',
              delta_str(current['Isp_vac'], stoich['Isp_vac']))
    st.metric('Isp at sea level (1.013 bar)  [s]',
              f'{current["Isp_sl"]:.1f}',
              delta_str(current['Isp_sl'], stoich['Isp_sl']))

with c3:
    st.metric('Characteristic velocity c*  [m/s]',
              f'{current["c_star"]:.0f}',
              delta_str(current['c_star'], stoich['c_star'], '{:+.0f}'))
    st.metric('Mean molecular weight  [kg/kmol]',
              f'{current["M_mol"]:.2f}',
              delta_str(current['M_mol'], stoich['M_mol'], '{:+.2f}'))

st.caption(f'Deltas are relative to the stoichiometric point '
           f'(O/F = {pa.stoichiometric_OF(prop["fuel"], prop["oxidizer"]):.3f}) '
           f'at the same chamber and exit pressure.')


# ---------------------------------------------------------------------------
# Curves
# ---------------------------------------------------------------------------

st.markdown('---')
col_a, col_b = st.columns(2)

with col_a:
    st.subheader(f'Isp vs O/F  —  {prop_name}')
    curve = of_sweep(prop_name, p_chamber, p_exit, n=30)
    fig_line = px.line(curve, x='OF_ratio', y='Isp_vac',
                       labels={'OF_ratio': 'O/F ratio [-]',
                               'Isp_vac': f'Isp at {p_exit:g} bar [s]'})
    fig_line.add_vline(x=OF, line_dash='dash', line_color='red',
                       annotation_text=f'current  OF={OF:.2f}',
                       annotation_position='top')
    of_stoich = pa.stoichiometric_OF(prop['fuel'], prop['oxidizer'])
    if of_min <= of_stoich <= of_max:
        fig_line.add_vline(x=of_stoich, line_dash='dot', line_color='gray',
                           annotation_text=f'stoich  OF={of_stoich:.2f}',
                           annotation_position='bottom right',
                           annotation_xshift=6)
    fig_line.update_layout(height=380, margin=dict(l=10, r=10, t=30, b=10))
    st.plotly_chart(fig_line, use_container_width=True)

with col_b:
    st.subheader('Chamber composition (top species)')
    sp_names = [current[f'species_{i}_name'] for i in range(1, pa.N_TOP_SPECIES + 1)]
    sp_fracs = [current[f'species_{i}_frac'] for i in range(1, pa.N_TOP_SPECIES + 1)]
    pairs = [(n, f) for n, f in zip(sp_names, sp_fracs) if n and f >= 0.01]
    if pairs:
        df_pie = pd.DataFrame(pairs, columns=['species', 'mole_fraction'])
        fig_pie = px.pie(df_pie, values='mole_fraction', names='species', hole=0.35)
        # Hide on-slice labels for species < 4 % — they still appear in the
        # legend on the right, but small wedges no longer get clipped text.
        slice_text = [f'{n}<br>{f * 100:.1f}%' if f >= 0.04 else ''
                      for n, f in pairs]
        fig_pie.update_traces(
            text=slice_text,
            textinfo='text',
            textposition='inside',
            insidetextorientation='radial',
            hovertemplate='<b>%{label}</b><br>mole fraction = %{value:.3f}<extra></extra>',
        )
        fig_pie.update_layout(height=380, margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(fig_pie, use_container_width=True)
    else:
        st.info('No species exceeds the 1% mole fraction threshold.')


# ---------------------------------------------------------------------------
# Real-engine validation table
# ---------------------------------------------------------------------------

st.markdown('---')
st.subheader('Validation — how the model compares to real engines')
st.caption('Calculation uses the same (O/F, p_chamber) as the real engine and '
           'expands to vacuum (p_exit = 0.001 bar). Real engines run below '
           'ideal Isp due to finite kinetics, boundary-layer losses, and finite '
           'nozzle expansion — a positive "Diff %" of 5–20% is expected.')
val = real_engine_table(p_exit_bar=0.001)
st.dataframe(val, use_container_width=True, hide_index=True)


# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------

with st.expander('About this dashboard'):
    st.markdown(
        '''
        - **Engine**: `propellant_analysis.py` — Cantera HP equilibrium for the
          chamber, isentropic + SP equilibrium for the nozzle exit.
        - **Mechanism**: GRI-Mech 3.0 (53 species).  Species substitutions:
          C2H5OH → C2H6, N2H4 → NH3, RP-1 → C2H6 (validation only).
        - **Caching**: per-propellant O/F sweeps are cached; sliders re-render
          instantly after the first computation.
        - **Limitations**: gri30 nitrogen chemistry is incomplete; N2O4/NH3
          results are illustrative, not predictive.  Real engines also include
          finite-rate kinetics, film cooling, and boundary-layer losses that
          this idealised equilibrium model neglects.
        '''
    )
