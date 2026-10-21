"""Run data processing, fitting, simulations and figure generation in order."""
from pathlib import Path
import json
import platform
import subprocess
import sys
from importlib.metadata import version
ROOT = Path(__file__).resolve().parent
# Use the active Python environment for every pipeline stage.
for script in ['literature_productivity.py','prepare_hplc.py','fit_pk.py','simulate.py','plot_results.py']:
    subprocess.run([sys.executable,str(ROOT/'src'/script)],check=True,cwd=ROOT)
metadata = {'python':platform.python_version(),
            'packages':{name:version(name) for name in ['numpy','scipy','matplotlib','scikit-learn']}}
(ROOT/'results/environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
print('Complete. Results are in results/.')
