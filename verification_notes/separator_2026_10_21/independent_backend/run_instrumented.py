from pathlib import Path
import json,subprocess,concurrent.futures,time,hashlib,re
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'xy_global_separator'/'distribution'
cases=json.loads((ROOT/'cases.json').read_text());out=ROOT/'results';out.mkdir(exist_ok=True)
def run_one(c):
 res=out/c['stem'];res.mkdir(exist_ok=True)
 with (res/'chart.json').open('w') as o,(res/'chart.log').open('w') as e:
  p=subprocess.run([str(ROOT/'chart/verify_guarded'),'1'],cwd=ROOT/'cases'/c['stem'],stdout=o,stderr=e)
 if p.returncode:raise RuntimeError((c['stem'],p.returncode,(res/'chart.log').read_text()))
 data=json.loads((res/'chart.json').read_text())
 if not(data['interior_pass'] and data['mu']>0):raise RuntimeError('guarded chart prerequisites failed')
 shutil_source=ROOT/'cases'/c['stem']/'results/coefficients.txt'
 (res/'coefficients.txt').write_bytes(shutil_source.read_bytes())
 exact=ROOT/'cases'/c['stem']/'results/coefficients_exact.hex'
 if not exact.is_file():raise RuntimeError('missing chart exact sidecar')
 (res/'coefficients_exact.hex').write_bytes(exact.read_bytes())
 return {'id':c['id'],'stem':c['stem'],'returncode':p.returncode,'guard_counters':(res/'chart.log').read_text().strip()}
t=time.time();results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for result in pool.map(run_one,cases):
  results.append(result)
  if len(results)%10==0:print('guarded charts:',len(results),flush=True)
with (out/'geometry.json').open('w') as o,(out/'geometry.log').open('w') as e:
 p=subprocess.run([str(ROOT/'geometry/verify_guarded'),str(SRC/'common/q.raw'),str(SRC/'geometry/cfg'),str(SRC/'centers')],cwd=ROOT/'geometry',stdout=o,stderr=e)
if p.returncode!=0:raise RuntimeError('guarded geometry failed')
exact=ROOT/'geometry/geometry_exact.hex'
if not exact.is_file():raise RuntimeError('missing exact geometry sidecar')
(out/'geometry_exact.hex').write_bytes(exact.read_bytes())
report={'scope':'diagnostic fixed-data replay with finite interval invariants; original arithmetic algorithm retained','charts':results,'geometry_returncode':p.returncode,'geometry_guard_counters':(out/'geometry.log').read_text().strip(),'duration_seconds':time.time()-t}
(out/'instrumented_report.json').write_text(json.dumps(report,indent=2)+'\n')
print('GUARDED_CHARTS_AND_GEOMETRY_PASS',flush=True)
