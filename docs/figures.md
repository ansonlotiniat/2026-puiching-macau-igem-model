# Figure and table index

| Output | Generator | Input |
| --- | --- | --- |
| `results/fitting/oral_fit.*` | `src/fit_pk.py` | Published serum means and SDs |
| `results/fitting/fit_results.json` | `src/fit_pk.py` | OLS/WLS optimization |
| `results/figures/hplc_productivity.*` | `src/plot_results.py` | Paired HPLC/OD rates |
| `results/figures/daily_dosing.*` | `src/plot_results.py` | 1/2/4 arrivals per day |
| `results/figures/OLS_1013_continuous_source.*` | `src/plot_results.py` | OLS source sweep |
| `results/figures/WLS_1013_continuous_source.*` | `src/plot_results.py` | WLS source sweep |
| `results/model/scenarios.csv` | `src/simulate.py` | Literature/culture productivity and configured parameters |

Each figure is exported as PNG, SVG and PDF. Numerical curve data are saved alongside the analysis results.
