# Methods

## HPLC and OD productivity

The `927` worksheet supplies peak area at 0 h (column D) and 24 h (column J). The calibration is:

```text
B [mM] = (peak_area - 10465) / 59906
Delta_B = B24 - B0
OD_integral [OD h] = 24 * (OD0 + OD24) / 2
qOD [mM / (OD h)] = Delta_B / OD_integral
q_v [mM/h per (CFU/mL)] = qOD / 8e8
```

OD is linearly interpolated between the two readings. Each culture condition has one paired observation; the OD dilution multiplier is 1. The per-cell conversion uses 8 × 10⁸ CFU/mL per OD unit. The resulting rate measures net accumulation in culture. Control contrasts subtract the separately OD-normalized control rate at the same initial butyrate concentration.

## Oral-dose pharmacokinetic fitting

Time is in minutes, concentrations in µg/mL and volumes in mL.

```text
dCg/dt = J/Vg - (kg + ka) Cg
dCb/dt = ka Cg Vg/Vb - ke Cb
```

The fitting input is six published sodium-butyrate group means. Set `J = 0`, `Vg = 1000`, `Vb = 5000` and `Cb(0) = 0.47`. Fit `[C0, kg, ka, ke]` from initial values `[20, 0.05, 0.05, 0.05]` with bounds `[(0,100), (0,0.5), (0,0.2), (0,0.1)]`.

`solve_ivp` uses RK45 with `rtol=1e-6`, `atol=1e-8`. L-BFGS-B minimizes mean squared error for OLS and the sum of squared residuals weighted by `1/SD²` for WLS. The JSON output records full-precision parameters, optimizer status, objective, RMSE and R². Weighted R² uses the weighted observation mean in its denominator.

## Retention and continuous production

The baseline daily input is 10¹⁰ CFU/day. Each arrival adds `F * daily_CFU / splits` cells, followed by exponential loss at 0.14/day. The model uses a 400 mL production volume and 1/hour luminal butyrate removal.

```text
mean_retained_CFU = F * daily_CFU / kN
production_mM_h = q_v * mean_retained_CFU / lumen_mL
lumen_increment_mM = production_mM_h / kB
J_ug_min = production_mM_h * lumen_mL * 88.11 / 60
serum_gain = ka / [Vb * ke * (kg + ka)]
serum_increment_ss = J_ug_min * serum_gain
```

For constant source activation, let `a=kg+ka` and `e=ke`. The serum increment at time `t` is its steady-state value multiplied by:

```text
1 - [a exp(-e t) - e exp(-a t)] / (a - e)
```

For each positive productivity, the daily input required for a target increment follows from linear scaling: `required_input = baseline_input * target / baseline_increment`. The configured comparison increments are 2.6 mM in the lumen and 0.10 µg/mL in serum. Nonpositive rates retain their signed forward values and have an empty inverse-dose field.

The parameter labels `OLS_1013`, `WLS_1013` and `OLS_legacy_wiki` identify the three report scenarios included in the HPLC calculation package. Each output row records its selected label.

## Literature productivity

`literature_productivity.py` converts the Bai 0.50 g/L endpoint using 88.1 g/mol, 48 hours, OD600=3.5 and 8 × 10⁸ CFU/mL/OD. The Wang input series is stored in `wang_culture.csv`; OD420 is divided by 1.5 and multiplied by the same CFU conversion. Trapezoidal integration supplies cumulative cell-hours. The 0–24 h concentration increase divided by cell-hours gives the growth productivity; the full-course slope is fitted through the origin.
