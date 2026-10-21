"""Plot HPLC productivity, retained cells and continuous-source serum responses."""
from pathlib import Path
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/figures'
OUT.mkdir(parents=True,exist_ok=True)
config = json.loads((ROOT/'config/parameters.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})

def read(path):
    with path.open(encoding='utf-8-sig') as handle: return list(csv.DictReader(handle))

def save(fig, name):
    for extension in ['png','svg','pdf']: fig.savefig(OUT/f'{name}.{extension}',dpi=220)
    plt.close(fig)

rows = read(ROOT/'results/hplc/paired_rates.csv')
fig, ax = plt.subplots(figsize=(7,4.5),layout='constrained')
for group,color in [('3+5','#8f2038'),('3+7','#347f75'),('Control','#777777')]:
    series = [r for r in rows if r['group']==group]
    # Connected points show the four measured culture conditions.
    ax.plot([float(r['added_mM']) for r in series],
            [float(r['qOD'])/config['CFU_per_mL_per_OD'] for r in series],'-o',color=color,label=group)
ax.axhline(0,color='#bbbbbb',linewidth=.8)
ax.set(xlabel='Added butyrate (mM)',ylabel='Apparent qᵥ [mM/h per (CFU/mL)]',title='HPLC and paired OD productivity')
ax.legend();save(fig,'hplc_productivity')
rows = read(ROOT/'results/model/daily_dosing.csv')
fig, ax = plt.subplots(figsize=(7,4.5),layout='constrained')
for n in [1,2,4]:
    ax.plot([float(r['time_day']) for r in rows],[float(r[f'{n}_doses_CFU_mL']) for r in rows],label=f'{n} doses/day')
ax.set(xlabel='Time (days)',ylabel='Retained cells (CFU/mL)',title='Equal daily input, split dosing');ax.legend();save(fig,'daily_dosing')
for label in ['OLS_1013','WLS_1013']:
    rows = read(ROOT/f'results/model/{label}_source_curves.csv')
    fig, ax = plt.subplots(figsize=(7,4.5),layout='constrained')
    for j in [.2,.4,.6,.8,1.]:
        ax.plot([float(r['time_min']) for r in rows],[float(r[f'J_{j:.1f}_ug_min']) for r in rows],label=f'{j:.1f} µg/min')
    ax.set(xlabel='Time (min)',ylabel='Serum increment (µg/mL)',title=f'Continuous production · {label}')
    ax.legend();save(fig,label+'_continuous_source')
print('Plots: PNG, SVG and PDF files written to results/figures.')
