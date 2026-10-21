"""Convert paired HPLC peak areas and OD readings into apparent productivity."""
from pathlib import Path
import csv
import math
from workbook import sheet_cells

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/hplc'
OUT.mkdir(parents=True, exist_ok=True)
cells = sheet_cells(ROOT / 'data/raw/hplc_927.xlsx', '927')
with (ROOT / 'data/raw/od_readings.csv').open(encoding='utf-8-sig') as handle:
    readings = list(csv.DictReader(handle))
# Each OD pair maps to the first HPLC record in its condition block.
pairs = {}
for reading in readings:
    key = (reading['strain'], float(reading['added_butyrate_mM']))
    pairs.setdefault(key, {})[int(reading['time_h'])] = reading
rows = []
for (group, added), times in pairs.items():
    row = int(times[0]['hplc_row'])
    area0, area24 = cells[f'D{row}'], cells[f'J{row}']
    # Workbook calibration: concentration in mM = (peak area - intercept) / slope.
    b0, b24 = (area0 - 10465) / 59906, (area24 - 10465) / 59906
    assert math.isclose(b0, cells[f'E{row}'], abs_tol=1e-10)
    assert math.isclose(b24, cells[f'K{row}'], abs_tol=1e-10)
    od0, od24 = float(times[0]['OD600']), float(times[24]['OD600'])
    # Linear interpolation between endpoints gives the 24-hour OD integral.
    od_area = 24 * (od0 + od24) / 2
    rows.append(dict(group=group, added_mM=added, row=row, delta_mM=b24-b0,
                     OD0=od0, OD24=od24, OD_area_linear=od_area,
                     qOD=(b24-b0)/od_area, paired_observations=1))
controls = {r['added_mM']: r['qOD'] for r in rows if r['group'] == 'Control'}
for row in rows:
    # Subtract productivity after normalizing each culture by its own OD integral.
    row['control_adjusted_qOD'] = row['qOD'] - controls[row['added_mM']]
with (OUT/'paired_rates.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
# Verify the reconstruction against the supplied analysis table.
with (ROOT/'data/processed/paired_rates.csv').open(encoding='utf-8-sig') as handle:
    expected = list(csv.DictReader(handle))
assert len(expected) == len(rows) == 12
for actual, reference in zip(rows, expected):
    assert actual['group'] == reference['group']
    for key in actual:
        if key != 'group':
            assert math.isclose(float(actual[key]),float(reference[key]),rel_tol=1e-11,abs_tol=1e-12), key
print('HPLC: 12 paired records reproduced from workbook peak areas and OD readings.')
