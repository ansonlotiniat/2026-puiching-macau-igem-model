# PuiChing Macau — Modeling code and data

Reproduce the HPLC/OD productivity analysis, oral sodium-butyrate fitting, daily cell-retention simulations and gut-to-blood calculations.

## Run the analysis

Use Python 3.14.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py
```

On Windows, activate the environment with `.venv\Scripts\activate`.

The pipeline writes CSV tables, parameter estimates and PNG/SVG/PDF figures to `results/`.

## Repository contents

| Path | Contents |
| --- | --- |
| `data/raw/` | HPLC workbook, OD source images and transcribed OD readings |
| `data/processed/paired_rates.csv` | Reference table of 12 paired culture conditions |
| `data/literature/` | Published serum means/SDs and literature productivity inputs |
| `config/parameters.json` | Volumes, arrival fractions, loss rates, source conversion and report PK parameter sets |
| `src/` | Data processing, fitting, simulation and plotting scripts |
| `notebooks/pk_fitting.ipynb` | Notebook entry point for the OLS/WLS fitting script |
| `docs/` | Methods, data dictionary and figure index |
| `results/` | Generated tables, figures and environment record |

## Analysis stages

`literature_productivity.py` first derives the Bai and Wang productivity inputs and exports the Wang cell-hour table.

1. `prepare_hplc.py` converts HPLC peak areas to concentration, pairs the 0/24 h OD readings and calculates apparent productivity.
2. `fit_pk.py` fits the six serum means using OLS and inverse-SD² WLS, preserving the supplied notebook's fourth-cell equations and optimizer settings.
3. `simulate.py` evaluates the named parameter sets in `config/parameters.json`, produces 132 scenarios and calculates continuous-source and split-dose trajectories.
4. `plot_results.py` exports publication figures in three formats.

Fitted estimates are saved in `results/fitting/fit_results.json`. The simulation stage selects its transport parameters from `config/parameters.json`; edit that file to run a different parameter scenario.

## Read the results

- [Methods and parameter definitions](docs/methods.md)
- [Input and output data dictionary](docs/data.md)
- [Figure and table index](docs/figures.md)
- [Source references](docs/references.md)

The GitHub Actions workflow reruns the pipeline for each push and saves the generated results as a downloadable artifact.

Project: [PuiChing Macau 2026 Model page](https://2026.igem.wiki/puiching-macau/model).
