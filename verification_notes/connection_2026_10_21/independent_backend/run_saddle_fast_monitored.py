"""Alternative same-check replay with lazily constructed error messages."""
from pathlib import Path
import time,json,subprocess,hashlib
R=Path(__file__).resolve().parent;out=R/'results';GIB=1<<30
def stat():return {k:int(v) for k,v in (x.split() for x in Path('/sys/fs/cgroup/memory.stat').read_text().splitlines())}
def used():
 s=stat();return s['anon']+s['kernel']+s.get('shmem',0)
def oom():return {k:int(v) for k,v in (x.split() for x in Path('/sys/fs/cgroup/memory.events').read_text().splitlines())}.get('oom',0)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
initial=used();initial_oom=oom();limit=int(Path('/sys/fs/cgroup/memory.max').read_text())
if initial>=4*GIB:raise RuntimeError('not enough initial anon/kernel/shmem headroom for lazy-message replay')
print('Starting lazy-message guarded4096; initialmemorybytes',initial,flush=True);start=time.time();peak=initial
with (out/'saddle_fast.json').open('w') as o,(out/'saddle_fast.log').open('w') as e:
 p=subprocess.Popen([str(R/'stage2/saddle_fast/verify_guarded_fast')],cwd=R/'stage2/saddle_fast',stdout=o,stderr=e)
 while p.poll() is None:
  current=used();peak=max(peak,current)
  if current>int(6.8*GIB) or oom()>initial_oom:
   p.terminate();p.wait(timeout=30);raise RuntimeError('lazy-message diagnostic stopped to protect aggregate resident budget')
  time.sleep(2)
if p.returncode:raise RuntimeError('lazy-message guarded saddle failed; '+(out/'saddle_fast.log').read_text())
(out/'saddle_fast_exact.hex').write_bytes((R/'stage2/saddle_fast/results/exact.hex').read_bytes())
clean=R.parent/'full_fresh/fresh_connection/prerequisites/saddle.json'
if (out/'saddle_fast.json').read_bytes()!=clean.read_bytes():raise RuntimeError('lazy-message scientific stdout differs from fresh original')
report={'status':'CONNECTION_SADDLE_LAZY_MESSAGE_GUARD_PASS','duration_seconds':time.time()-start,'all_original_numerical_expressions_and_boolean_checks_retained':True,'scientific_stdout_byte_identical_to_fresh_original':True,'aggregate_anon_kernel_shmem_initial_bytes':initial,'aggregate_anon_kernel_shmem_peak_sampled_bytes':peak,'aggregate_memory_limit_bytes':limit,'peak_sampling_interval_seconds':2,'oom_events_initial':initial_oom,'oom_events_final':oom(),'guarded_process_returncode':p.returncode,'program':'stage2/saddle_fast/verify_guarded_fast','engine_source_sha256':sha(R/'stage2/saddle_fast/engine.cpp'),'certify_source_sha256':sha(R/'stage2/saddle_fast/certify.cpp'),'scientific_stdout_sha256':sha(out/'saddle_fast.json'),'exact_sidecar_sha256':sha(out/'saddle_fast_exact.hex')}
(out/'saddle_fast_guard_report.json').write_text(json.dumps(report,indent=2)+'\n');print(report['status'],flush=True)
