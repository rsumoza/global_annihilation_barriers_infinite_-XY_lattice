"""Run the bounded scalar probe from exact completed baseline sidecar values."""
from pathlib import Path
import subprocess,json
R=Path(__file__).resolve().parent;prior=R.parent.parent/'audit_separator_2026_10_21/independent_backend/results/baseline.log'
line=next(x for x in prior.read_text().splitlines() if x.startswith('EXACT_BASELINE '));v=dict(w.split('=',1) for w in line.split()[1:])
flags=json.loads((R/'instrumentation.json').read_text())[0]['flags']
subprocess.run(['g++',*flags,'baseline_scalar_probe.cpp','-o','baseline_scalar_probe'],cwd=R,check=True)
with (R/'results/baseline_scalar_exact.hex').open('w') as o,(R/'results/baseline_scalar_probe.log').open('w') as e:
 p=subprocess.run([str(R/'baseline_scalar_probe'),v['c'],v['rho']],stdout=o,stderr=e)
if p.returncode:raise RuntimeError('baseline scalar probe failed')
print('BASELINE_SCALAR_PROBE_BOUND_TO_PRIOR_EXACT_VALUES')
