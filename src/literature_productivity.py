"""Calculate literature productivity from an endpoint and a culture time course."""
from pathlib import Path
import json
import math
import numpy as np
from scipy.integrate import cumulative_trapezoid
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/literature'
OUT.mkdir(parents=True,exist_ok=True)
data = np.genfromtxt(ROOT/'data/literature/wang_culture.csv',delimiter=',',names=True)
# Convert optical density at 420 nm to OD600 before converting to cell density.
cells = data['OD420']/1.5*8e8
integral = cumulative_trapezoid(cells,data['time_h'],initial=0)
change = data['butyrate_mM']-data['butyrate_mM'][0]
# Growth productivity uses the concentration increase over the 0–24 h integral.
index = int(np.flatnonzero(data['time_h']==24)[0])
wang = float(change[index]/integral[index])
# The full-course regression is constrained to pass through the origin.
slope = float(np.dot(integral,change)/np.dot(integral,integral))
bai = (.5/88.1*1000)/(48*3.5*8e8)
result = {'Bai_q_v':bai,'Wang_growth_q_v':wang,'Wang_full_course_slope':slope}
reference = json.loads((ROOT/'data/literature/part_i_source.json').read_text())['q_v_table']
assert math.isclose(bai,reference['Bai_Mansell_2020_EcN']['glucose_central_anchor']['q_v_central'],rel_tol=1e-12)
assert math.isclose(wang,reference['Wang_2019_Ecoli_W']['q_growth_upper_bound'],rel_tol=1e-12)
(OUT/'productivity.json').write_text(json.dumps(result,indent=2)+'\n')
np.savetxt(OUT/'cell_hours.csv',np.column_stack([data['time_h'],cells,integral,change]),delimiter=',',
           header='time_h,CFU_mL,cell_hours_CFU_h_mL,butyrate_change_mM',comments='')
print('Literature: Bai endpoint and Wang cell-hour productivity reproduced.')
