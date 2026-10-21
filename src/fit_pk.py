"""Fit oral sodium-butyrate group means with OLS and WLS (notebook cell 4)."""
from pathlib import Path
import json
import csv
import matplotlib
matplotlib.use("Agg")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/fitting"
OUT.mkdir(parents=True, exist_ok=True)
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from sklearn.metrics import r2_score
records = np.genfromtxt(ROOT / 'data/literature/serum_nab.csv', delimiter=',', names=True)
t_data = records['time_min']
c_blood_observed = records['mean_ug_mL']
std_observed = records['sd_ug_mL']

# Concentrations use micrograms/mL and time uses minutes.
V_gut = 1000.0
V_blood = 5000.0
V_ratio_1 = V_gut / V_blood

# Transfer conserves amount through the gut-to-blood volume ratio.
def gbb_ode_system(t, y, params):
    C_gut, C_blood = y
    v_produced, v_gutremoved, v_absorption, v_bloodremoved = params
    dC_gut = v_produced / V_gut - (v_gutremoved + v_absorption) * C_gut
    dC_blood = v_absorption * C_gut * V_ratio_1 - v_bloodremoved * C_blood
    return [dC_gut, dC_blood]

# OLS assigns equal weight to each of the six concentration means.
def loss_mse(free_params):
    C_gut_init, v_gutremoved, v_absorption, v_bloodremoved = free_params
    if any((p < 0 for p in free_params)):
        return 1000000.0
    y0 = [C_gut_init, c_blood_observed[0]]
    sol = solve_ivp(gbb_ode_system, [t_data[0], t_data[-1]], y0, args=([0.0, v_gutremoved, v_absorption, v_bloodremoved],), t_eval=t_data, method='RK45', rtol=1e-06, atol=1e-08)
    if sol.status != 0 or len(sol.t) < len(t_data):
        return 1000000.0
    return np.mean((sol.y[1] - c_blood_observed) ** 2)

# WLS uses inverse squared published SD as in notebook cell 4.
def loss_wls(free_params):
    C_gut_init, v_gutremoved, v_absorption, v_bloodremoved = free_params
    if any((p < 0 for p in free_params)):
        return 1000000.0
    y0 = [C_gut_init, c_blood_observed[0]]
    sol = solve_ivp(gbb_ode_system, [t_data[0], t_data[-1]], y0, args=([0.0, v_gutremoved, v_absorption, v_bloodremoved],), t_eval=t_data, method='RK45', rtol=1e-06, atol=1e-08)
    if sol.status != 0 or len(sol.t) < len(t_data):
        return 1000000.0
    weights = 1.0 / std_observed ** 2
    return np.sum(weights * (sol.y[1] - c_blood_observed) ** 2)
# Parameter order: initial gut concentration, gut loss, absorption, blood loss.
initial_guess = [20.0, 0.05, 0.05, 0.05]
bounds_intermediate = [(0, 100), (0, 0.5), (0, 0.2), (0, 0.1)]
# Run bounded L-BFGS-B with the original SciPy default stopping settings.
res_mse = minimize(loss_mse, initial_guess, bounds=bounds_intermediate, method='L-BFGS-B')
res_wls = minimize(loss_wls, initial_guess, bounds=bounds_intermediate, method='L-BFGS-B')
p_mse = res_mse.x
sol_mse = solve_ivp(gbb_ode_system, [t_data[0], t_data[-1]], [p_mse[0], c_blood_observed[0]], args=([0.0, p_mse[1], p_mse[2], p_mse[3]],), t_eval=t_data, rtol=1e-06, atol=1e-08)
fitted_mse = sol_mse.y[1]
rmse_mse = np.sqrt(np.mean((fitted_mse - c_blood_observed) ** 2))
r2_unweighted_mse = r2_score(c_blood_observed, fitted_mse)
p_wls = res_wls.x
sol_wls = solve_ivp(gbb_ode_system, [t_data[0], t_data[-1]], [p_wls[0], c_blood_observed[0]], args=([0.0, p_wls[1], p_wls[2], p_wls[3]],), t_eval=t_data, rtol=1e-06, atol=1e-08)
fitted_wls = sol_wls.y[1]
rmse_wls = np.sqrt(np.mean((fitted_wls - c_blood_observed) ** 2))
r2_unweighted_wls = r2_score(c_blood_observed, fitted_wls)
weights = 1.0 / std_observed ** 2
y_weighted_mean = np.sum(weights * c_blood_observed) / np.sum(weights)
ss_res_w = np.sum(weights * (c_blood_observed - fitted_wls) ** 2)
ss_tot_w = np.sum(weights * (c_blood_observed - y_weighted_mean) ** 2)
r2_weighted_wls = 1.0 - ss_res_w / ss_tot_w


# Save full-precision estimates, convergence information and observation residuals.
results = {}
for name, result, prediction in [('OLS', res_mse, fitted_mse), ('WLS', res_wls, fitted_wls)]:
    if not result.success:
        raise RuntimeError(f"{name} optimization failed: {result.message}")
    results[name] = {'parameters': dict(zip(['C0', 'kg', 'ka', 'ke'], map(float, result.x))),
                     'objective': float(result.fun), 'success': bool(result.success),
                     'message': str(result.message),
                     'RMSE': float(np.sqrt(np.mean((prediction-c_blood_observed)**2))),
                     'R2_unweighted': float(r2_score(c_blood_observed, prediction))}
results['WLS']['R2_weighted'] = float(r2_weighted_wls)
(OUT / 'fit_results.json').write_text(json.dumps(results, indent=2) + '\n')
with (OUT / 'predictions.csv').open('w', newline='') as handle:
    writer = csv.writer(handle)
    writer.writerow(['time_min','observed_ug_mL','SD_ug_mL','OLS_ug_mL','WLS_ug_mL','OLS_residual','WLS_residual'])
    writer.writerows(zip(t_data,c_blood_observed,std_observed,fitted_mse,fitted_wls,
                         fitted_mse-c_blood_observed,fitted_wls-c_blood_observed))
# Evaluate the same ODE on a dense time grid for plotting.
time = np.linspace(0, 210, 421)
fig, ax = plt.subplots(figsize=(8, 5), layout='constrained')
ax.errorbar(t_data,c_blood_observed,yerr=std_observed,fmt='o',color='#343434',capsize=3,label='Mean ± SD')
curves = []
for label, result, color in [('OLS',res_mse,'#8f2038'),('WLS',res_wls,'#347f75')]:
    p = result.x
    sol = solve_ivp(gbb_ode_system,[0,210],[p[0],c_blood_observed[0]],
                    args=([0.,p[1],p[2],p[3]],),t_eval=time,rtol=1e-6,atol=1e-8)
    if not sol.success: raise RuntimeError(sol.message)
    ax.plot(time,sol.y[1],color=color,label=label)
    curves.append(sol.y[1])
ax.set(xlabel='Time (min)',ylabel='Serum butyrate (µg/mL)',title='Oral sodium-butyrate fit')
ax.legend(); ax.spines[['top','right']].set_visible(False)
for ext in ['png','svg','pdf']: fig.savefig(OUT / f'oral_fit.{ext}',dpi=220)
plt.close(fig)
np.savetxt(OUT/'curves.csv',np.column_stack([time,*curves]),delimiter=',',
           header='time_min,OLS_ug_mL,WLS_ug_mL',comments='')
print('Fitting: OLS and WLS estimates, residuals and figures written to results/fitting.')
