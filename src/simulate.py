"""Calculate retention, production, serum responses and inverse daily inputs."""
from pathlib import Path
import csv
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/model'
OUT.mkdir(parents=True, exist_ok=True)
config = json.loads((ROOT/'config/parameters.json').read_text())

def gain(pk):
    """Steady serum increment per unit continuous source [min/mL]."""
    return pk['ka'] / (config['blood_mL'] * pk['ke'] * (pk['kg'] + pk['ka']))

def forward(dose, q, arrival, pk):
    """Map CFU/day and cellular productivity to lumen and serum increments."""
    retained = arrival * dose / config['cell_loss_per_day']
    production = q * retained / config['lumen_mL']
    # mM * mL equals micromoles; molecular mass converts to micrograms.
    source = production * config['lumen_mL'] * config['butyrate_molar_mass_g_mol'] / 60
    return retained, production, production/config['gut_removal_per_h'], source, source*gain(pk)

def step(time, source, pk):
    """Analytical serum response to a constant source switched on at time zero."""
    a, e = pk['kg'] + pk['ka'], pk['ke']
    fraction = 1 - (a*np.exp(-e*time) - e*np.exp(-a*time))/(a-e)
    return source * gain(pk) * fraction

def write_csv(name, rows):
    with (OUT/name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

qtable = json.loads((ROOT/'data/literature/part_i_source.json').read_text())['q_v_table']
scenarios = [dict(name='Bai',kind='literature',q=qtable['Bai_Mansell_2020_EcN']['glucose_central_anchor']['q_v_central']),
             dict(name='Wang',kind='literature',q=qtable['Wang_2019_Ecoli_W']['q_growth_upper_bound'])]
with (ROOT/'results/hplc/paired_rates.csv').open() as handle:
    for row in csv.DictReader(handle):
        name = f"{row['group']} / {float(row['added_mM']):g} mM"
        scenarios.append(dict(name=name,kind='apparent_net',q=float(row['qOD'])/config['CFU_per_mL_per_OD']))
        if row['group'] != 'Control':
            scenarios.append(dict(name=name,kind='control_contrast',q=float(row['control_adjusted_qOD'])/config['CFU_per_mL_per_OD']))
rows = []
for scenario in scenarios:
    q = scenario['q']
    for formulation, arrival in config['arrival_fractions'].items():
        for label, pk in config['pk_sets'].items():
            retained, production, lumen, source, serum = forward(config['daily_CFU'],q,arrival,pk)
            # Linear scaling makes inverse inputs proportional to the target increment.
            gut_dose = config['daily_CFU']*config['gut_target_mM']/lumen if q > 0 else None
            serum_dose = config['daily_CFU']*config['serum_target_ug_mL']/serum if q > 0 else None
            if q > 0:
                assert math.isclose(forward(gut_dose,q,arrival,pk)[2],config['gut_target_mM'],rel_tol=1e-12)
                assert math.isclose(forward(serum_dose,q,arrival,pk)[4],config['serum_target_ug_mL'],rel_tol=1e-12)
            rows.append(dict(scenario=scenario['name'],kind=scenario['kind'],q_v=q,F=arrival,fit=label,
                             daily_CFU=config['daily_CFU'],retained_CFU=retained,density_CFU_mL=retained/config['lumen_mL'],
                             production_mM_h=production,lumen_increment_mM=lumen,source_ug_min=source,
                             serum_increment_ss_ug_mL=serum,serum_increment_210min_constant_source_ug_mL=float(step(210,source,pk)),
                             gut_target_mM=config['gut_target_mM'],serum_target_ug_mL=config['serum_target_ug_mL'],
                             gut_required_CFU_day=gut_dose,serum_required_CFU_day=serum_dose))
write_csv('scenarios.csv', rows)
# Save the full 0–210 minute source sweep for both report parameter sets.
time = np.arange(0,210.5,.5)
for label in ['OLS_1013','WLS_1013']:
    sweep = [{'time_min':float(t),**{f'J_{j:.1f}_ug_min':float(step(t,j,config['pk_sets'][label])) for j in [.2,.4,.6,.8,1.]}} for t in time]
    write_csv(label+'_source_curves.csv',sweep)
# The daily dose is split into one, two or four instantaneous arrivals.
time_days = np.linspace(0,21,4201)
dosing = {'time_day':time_days}
for splits in [1,2,4]:
    population = np.zeros_like(time_days)
    for pulse in np.arange(0,21,1/splits):
        elapsed = time_days-pulse
        population += np.where(elapsed>=0,.5*config['daily_CFU']/splits*np.exp(-config['cell_loss_per_day']*np.maximum(elapsed,0)),0)
    dosing[f'{splits}_doses_CFU_mL'] = population/config['lumen_mL']
write_csv('daily_dosing.csv',[dict(zip(dosing,values)) for values in zip(*dosing.values())])
print(f'Model: {len(rows)} scenarios, source sweeps and daily-dosing curves written.')
