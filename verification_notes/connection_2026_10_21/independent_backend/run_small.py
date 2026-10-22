from pathlib import Path
import json,subprocess,time,concurrent.futures
R=Path(__file__).resolve().parent;S=R.parent/'xy_global_connection/distribution';out=R/'results';out.mkdir(exist_ok=True)
def run(name,cmd,cwd):
 with (out/(name+'.json')).open('w') as o,(out/(name+'.log')).open('w') as e:
  p=subprocess.run(list(map(str,cmd)),cwd=cwd,stdout=o,stderr=e)
 if p.returncode:raise RuntimeError('guarded stage failed '+name+'; '+(out/(name+'.log')).read_text())
 return {'program':name,'returncode':p.returncode}
def paths():
 ret=[]
 for name,flag in [('leg_return','--leg'),('leg_annihilation','--terminal-leg'),('full_annihilation',None)]:
  cmd=[R/'connections/path_guarded',R/'connections'/name,out/name]
  if flag:cmd.append(flag)
  ret.append(run(name,cmd,R/'connections'))
 return ret
start=time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 b=pool.submit(run,'barrier',[R/'stage2/barrier/verify_guarded','640'],R/'stage2/barrier')
 ps=pool.submit(paths)
 results=[b.result(),*ps.result()]
results.append(run('join_points',[R/'connections/join_guarded',S/'stage2/saddle/inputs/h.raw',S/'stage2/common/base_h.raw',S/'stage2/saddle/inputs/q.raw',R/'connections/c_minus.z.raw',R/'connections/c_plus.z.raw','256'],R/'connections'))
(out/'barrier_exact.hex').write_bytes((R/'stage2/barrier/results/exact.hex').read_bytes());(out/'join_exact.hex').write_bytes((R/'connections/results/join_exact.hex').read_bytes())
(out/'small_guard_report.json').write_text(json.dumps({'status':'CONNECTION_SMALL_GUARD_PASS','programs':results,'duration_seconds':time.time()-start},indent=2)+'\n')
print('CONNECTION_SMALL_GUARD_PASS',flush=True)
