# Data dictionary

## Inputs

| File | Contents and units |
| --- | --- |
| `data/raw/hplc_927.xlsx` | Worksheet `927`: sample IDs, peak areas and calibrated butyrate concentrations |
| `data/raw/od_0h.jpg`, `od_24h.jpg` | Source images for the paired OD600 measurements |
| `data/raw/od_readings.csv` | Strain, initial added butyrate (mM), time (h), OD600, dilution multiplier and HPLC row |
| `data/processed/paired_rates.csv` | Twelve conditions: 3+5, 3+7 and Control at 0/10/20/50 mM added butyrate |
| `data/literature/serum_nab.csv` | Time (min), serum mean and SD (µg/mL), participant count |
| `data/literature/part_i_source.json` | Literature productivity values and associated Part I calculations |

HPLC row mapping: 3+5 → 2/8/14/20; 3+7 → 27/33/39/45; Control → 51/57/63/69, in increasing added-butyrate order.

`paired_rates.csv` fields: `delta_mM` is the 24-hour concentration increase; `OD_area_linear` is the trapezoidal OD integral in OD·h; `qOD` is the increase divided by that integral; `control_adjusted_qOD` subtracts the matched control's qOD. `paired_observations` is 1 for every condition.

## Outputs

`results/model/scenarios.csv` contains 22 productivity scenarios × 2 arrival fractions × 3 transport parameter sets = 132 rows. `kind` distinguishes literature, apparent net productivity and normalized control contrast. Concentration outputs are identified by `mM` or `ug_mL`; `source_ug_min` is the amount input to the PK model. Required daily inputs are in CFU/day.

The source-sweep CSVs use a 0.5-minute grid from 0 to 210 minutes. Each `J_*_ug_min` column contains the predicted serum increment in µg/mL for the source named in the header.

`results/fitting/predictions.csv` contains the six observations, fitted values and signed residuals (prediction minus observation). `curves.csv` contains the dense plotting grid. `results/environment.json` records Python and package versions.
