"""Complete sidecars for a run begun before exact exporters were installed.
Not required for the documented fresh workflow (install exporters before running).
"""
from pathlib import Path
import subprocess,json,time,concurrent.futures
ROOT=Path(__file__).resolve().parent
deadline=time.monotonic()+2400
while not (ROOT/'results/instrumented_report.json').exists():
 if time.monotonic()>deadline:raise RuntimeError('guarded first pass did not finish')
 time.sleep(5)
cases=json.loads((ROOT/'cases.json').read_text());missing=[c for c in cases if not (ROOT/'cases'/c['stem']/'results/coefficients_exact.hex').exists()]
def one(c):
 out=ROOT/'results'/c['stem'];inp=ROOT/'cases'/c['stem']
 with (out/'chart.json').open('w') as o,(out/'chart.log').open('w') as e:
  p=subprocess.run([str(ROOT/'chart/verify_guarded'),'1'],cwd=inp,stdout=o,stderr=e)
 if p.returncode!=0:raise RuntimeError('sidecar replay failed '+c['stem'])
 (out/'coefficients.txt').write_bytes((inp/'results/coefficients.txt').read_bytes())
 return c['id']
print('Replaying missing exact-sidecar cases:',len(missing),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for k,result in enumerate(pool.map(one,missing),1):
  if k%10==0:print('sidecar replays',k,flush=True)
for c in cases:
 exact=ROOT/'cases'/c['stem']/'results/coefficients_exact.hex'
 if not exact.is_file():raise RuntimeError('missing exact case sidecar')
 (ROOT/'results'/c['stem']/'coefficients_exact.hex').write_bytes(exact.read_bytes())
exact=ROOT/'geometry/geometry_exact.hex'
if not exact.exists():
 historical=ROOT.parent.parent/'geometry_exact.hex'
 if not historical.is_file():raise RuntimeError('missing historical geometry sidecar')
 exact.write_bytes(historical.read_bytes())
(ROOT/'results/geometry_exact.hex').write_bytes(exact.read_bytes())
(ROOT/'results/sidecar_completion.json').write_text(json.dumps({'status':'SIDECAR_COMPLETION_PASS','replayed_case_ids':[c['id'] for c in missing]},indent=2)+'\n')
print('SIDECAR_COMPLETION_PASS',flush=True)
