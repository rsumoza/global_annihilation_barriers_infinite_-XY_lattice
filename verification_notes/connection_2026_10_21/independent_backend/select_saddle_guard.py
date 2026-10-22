"""Bind a completed guarded saddle to canonical diagnostic outputs."""
from pathlib import Path
import argparse,json,hashlib
R=Path(__file__).resolve().parent;out=R/'results'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
p=argparse.ArgumentParser();p.add_argument('--variant',choices=['original','lazy'],default='original');p.add_argument('--incomplete-extra-stopped',action='store_true');a=p.parse_args()
fast=a.variant=='lazy';stem='saddle_fast' if fast else 'saddle';source='stage2/saddle_fast' if fast else 'stage2/saddle';report=out/(stem+'_guard_report.json')
j=json.loads(report.read_text());expected='CONNECTION_SADDLE_LAZY_MESSAGE_GUARD_PASS' if fast else 'CONNECTION_SADDLE_GUARD_PASS';need(j['status']==expected,'selected guard run did not complete successfully')
clean=R.parent/'full_fresh/fresh_connection/prerequisites/saddle.json';need((out/(stem+'.json')).read_bytes()==clean.read_bytes(),'selected saddle stdout differs from fresh original')
if fast:
 for src,dst in [(stem+'.json','saddle.json'),(stem+'.log','saddle.log'),(stem+'_exact.hex','saddle_exact.hex')]:
  (out/dst).write_bytes((out/src).read_bytes())
selected={'status':'COMPLETED_SADDLE_GUARD_SELECTED','variant':'lazy_messages' if fast else 'original_eager_message_guard','source_directory':source,'program':source+('/verify_guarded_fast' if fast else '/verify_guarded'),'engine_source_sha256':sha(R/source/'engine.cpp'),'certify_source_sha256':sha(R/source/'certify.cpp'),'stdout_sha256':sha(out/'saddle.json'),'exact_sidecar_sha256':sha(out/'saddle_exact.hex'),'completed_guard_report_sha256':sha(report),'incomplete_alternative_stopped_and_excluded':a.incomplete_extra_stopped}
(R/'selected_saddle_guard_binding.json').write_text(json.dumps(selected,indent=2)+'\n');print('COMPLETED_SADDLE_GUARD_SELECTED',a.variant)
