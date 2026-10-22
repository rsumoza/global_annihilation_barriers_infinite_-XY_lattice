"""Fixed4096 diagnostic replay with explicit aggregate memory headroom checks."""
from pathlib import Path
import time,json,subprocess
R=Path(__file__).resolve().parent;out=R/'results';GIB=1<<30
mempath=Path('/sys/fs/cgroup/memory.current');maxpath=Path('/sys/fs/cgroup/memory.max')
def stat():return {k:int(v) for k,v in (x.split() for x in Path('/sys/fs/cgroup/memory.stat').read_text().splitlines())}
def used():
 s=stat();return s['anon']+s['kernel']+s.get('shmem',0)
def oom():return {k:int(v) for k,v in (x.split() for x in Path('/sys/fs/cgroup/memory.events').read_text().splitlines())}.get('oom',0)
limit=int(maxpath.read_text());initial=used();initial_oom=oom()
if initial>=4*GIB:raise RuntimeError('not enough anon/kernel/shmem memory headroom; guarded saddle not started')
print('Starting guarded4096; initialmemorybytes',initial,flush=True);start=time.time();peak=initial
with (out/'saddle.json').open('w') as o,(out/'saddle.log').open('w') as e:
 p=subprocess.Popen([str(R/'stage2/saddle/verify_guarded')],cwd=R/'stage2/saddle',stdout=o,stderr=e)
 while p.poll() is None:
  current=used();peak=max(peak,current)
  if current>int(6.8*GIB) or oom()>initial_oom:
   p.terminate();p.wait(timeout=30);raise RuntimeError('guarded saddle stopped to protect aggregate memory headroom')
  time.sleep(2)
if p.returncode:raise RuntimeError('guarded saddle failed; '+(out/'saddle.log').read_text())
(out/'saddle_exact.hex').write_bytes((R/'stage2/saddle/results/exact.hex').read_bytes())
(out/'saddle_guard_report.json').write_text(json.dumps({'status':'CONNECTION_SADDLE_GUARD_PASS','duration_seconds':time.time()-start,'aggregate_anon_kernel_shmem_initial_bytes':initial,'aggregate_anon_kernel_shmem_peak_sampled_bytes':peak,'aggregate_memory_limit_bytes':limit,'peak_sampling_interval_seconds':2,'guarded_process_returncode':p.returncode,'original_overlap_authorized_after_headroom_check':True,'memory_current_not_used_as_resident_budget':'disk page cache is reclaimable','oom_events_initial':initial_oom,'oom_events_final':oom()},indent=2)+'\n')
print('CONNECTION_SADDLE_GUARD_PASS',flush=True)
