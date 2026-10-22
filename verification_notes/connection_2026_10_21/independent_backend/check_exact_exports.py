"""Exact fixed-data serialization audit; no binary float parser is used.

This diagnostic retains the original arithmetic formulas and does not claim a
second interval arithmetic implementation. Only freshly rebuilt connection JSON
is accepted as the reference scientific output.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,re,sys
# Exact finite extended64 subnormals have rational denominators with about5000
# decimal digits; this budget affects report serialization, not arithmetic.
sys.set_int_max_str_digits(30000)
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'xy_global_connection/distribution'
CLEAN=ROOT.parent/'full_fresh/fresh_connection'
PAD=F('1e-18')
def need(ok,why):
 if not ok:raise RuntimeError(why)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(),parse_float=F)
def hf(s):
 sign=-1 if s.startswith('-') else 1;s=s.lstrip('+-');mant,e=s.lower().split('p');need(mant.startswith('0x'),'invalid hexfloat');a,sep,b=mant[2:].partition('.')
 q=F(sign*int(a+b,16),16**len(b));e=int(e)
 return q*(2**e if e>=0 else F(1,2**(-e)))
def at(obj,path):
 for k in path.split('/'):obj=obj[int(k)] if isinstance(obj,list) else obj[k]
 return obj

def one_bound(raw,printed,direction,label):
 raw=hf(raw);printed=F(printed);err=abs(raw-printed)
 need(err<PAD,label+': PAD insufficient')
 if direction=='L':need(printed-PAD<=raw,label+': lower direction failed')
 elif direction=='U':need(printed+PAD>=raw,label+': upper direction failed')
 else:raise RuntimeError('invalid direction')
 return {'field':label,'direction':direction,'raw_exact':str(raw),'printed_exact':str(printed),'absolute_error_exact':str(err),'pad_minus_error_exact':str(PAD-err)}

def summary(name,sidecar,cleanpath):
 p=ROOT/'results'/f'{name}.json';obj=load(p);need(p.read_bytes()==cleanpath.read_bytes(),name+': scientific stdout differs from fresh original')
 checks=[];values={};rawconstants={}
 for line in sidecar.read_text().splitlines():
  v=line.split()
  if v[0] in ('NODE','SEG'):continue
  if v[0]=='RAW':rawconstants[v[1]]=hf(v[2]);continue
  if v[0]=='I':
   need(len(v)==4,'malformed interval sidecar');need(hf(v[2])<=hf(v[3]),'unordered binary interval');q=at(obj,v[1]);need(len(q)==2,'interval shape')
   for k,dr in ((0,'L'),(1,'U')):checks.append(one_bound(v[k+2],q[k],dr,name+':'+v[1]+'/'+str(k)));values[v[1]+'/'+str(k)]=hf(v[k+2])
  else:
   need(len(v)==3,'malformed scalar sidecar');checks.append(one_bound(v[2],at(obj,v[1]),v[0],name+':'+v[1]));values[v[1]]=hf(v[2])
 return {'name':name,'summary_sha256':sha(p),'exact_sidecar_sha256':sha(sidecar),'checks':checks},values,rawconstants

def make_baseline_reuse():
 # Previously completed guarded separator baseline is identical source+input.
 prior=ROOT.parent.parent/'audit_separator_2026_10_21/independent_backend'
 text=(prior/'results/baseline.log').read_text();line=next(x for x in text.splitlines() if x.startswith('EXACT_BASELINE '));v=dict(w.split('=',1) for w in line.split()[1:])
 out=ROOT/'results';(out/'baseline.json').write_bytes((prior/'results/baseline.json').read_bytes());(out/'baseline_reused_guard.log').write_bytes((prior/'results/baseline.log').read_bytes())
 probe=(out/'baseline_scalar_exact.hex').read_text()
 (out/'baseline_exact.hex').write_text('U delta '+v['delta']+'\nU rho '+v['rho']+'\nL c '+v['c']+'\nL energy_candidate/0 '+v['E_lower']+'\n'+probe)
 # Verify every fixed minimum input, not merely h.raw or the source alone.
 binding=json.loads((ROOT/'baseline_reuse_binding.json').read_text())
 for row in binding['files']:need(sha(ROOT/'stage2/baseline'/row['file'])==row['sha256'],'reused minimum proposal changed')
 need(sha(SRC/'stage2/baseline/engine.cpp')==binding['engine_sha256'],'reused baseline original source changed')
 need(sha(out/'baseline_reused_guard.log')==binding['prior_guarded_log_sha256'],'prior guard record changed')
 return binding

def working_bindings():
 need(sha(SRC/'MANIFEST.sha256')=='ba100e9b7b8763bcd9ce1bcb7ee1a50a1f76669d2be92966ecbe6f6208feb27c','frozen connection manifest differs')
 rows=[]
 for line in (SRC/'MANIFEST.sha256').read_text().splitlines():
  digest,name=line.split('  ',1);rel=Path(name)
  need(not rel.is_absolute() and '..' not in rel.parts,'unsafe scientific manifest path')
  need(sha(SRC/rel)==digest,'frozen scientific member changed '+name)
  target=ROOT/rel
  if target.is_file() and target.suffix!='.cpp' and target.name!='MANIFEST.sha256':
   need(sha(target)==digest,'working fixed input differs '+name);rows.append({'path':name,'sha256':digest})
 base=sha(SRC/'stage2/common/base_h.raw')
 for name in ('stage2/baseline/min.h','stage2/saddle/inputs/reference_h.raw','stage2/barrier/inputs/h.raw','connections/leg_return.reference_h.raw','connections/leg_annihilation.reference_h.raw','connections/full_annihilation.reference_h.raw'):
  need(sha(ROOT/name)==base,'shared reference alias differs '+name);rows.append({'path':name,'sha256':base})
 n=256*256*8
 with (ROOT/'connections/leg_return.z.raw').open('rb') as f:f.seek(-n,2);plus=f.read(n)
 with (ROOT/'connections/leg_annihilation.z.raw').open('rb') as f:minus=f.read(n)
 for name,data in [('connections/c_plus.z.raw',plus),('connections/c_minus.z.raw',minus)]:
  need((ROOT/name).read_bytes()==data,'join endpoint not extracted from certified fixed leg '+name);rows.append({'path':name,'sha256':sha(ROOT/name)})
 fast=ROOT/'stage2/saddle_fast/inputs'
 if fast.is_dir():
  for p in (ROOT/'stage2/saddle/inputs').iterdir():
   if p.is_file():
    need(sha(fast/p.name)==sha(p),'lazy-message saddle input differs '+p.name);rows.append({'path':'stage2/saddle_fast/inputs/'+p.name,'sha256':sha(p)})
 return rows

def main():
 out=ROOT/'results';reuse=make_baseline_reuse();inputbindings=working_bindings()
 selected=load(ROOT/'selected_saddle_guard_binding.json')
 need(selected['status']=='COMPLETED_SADDLE_GUARD_SELECTED','no completed guarded saddle selected')
 need(sha(ROOT/selected['source_directory']/'engine.cpp')==selected['engine_source_sha256'],'selected saddle engine changed')
 need(sha(ROOT/selected['source_directory']/'certify.cpp')==selected['certify_source_sha256'],'selected saddle source changed')
 need(sha(out/'saddle.json')==selected['stdout_sha256'] and sha(out/'saddle_exact.hex')==selected['exact_sidecar_sha256'],'selected saddle output binding differs')
 need('GUARD_DIAGNOSTIC finite_interval_constructions=' in (out/'saddle.log').read_text(),'selected saddle guard log missing')
 cert=load(CLEAN/'certificate.json');need(cert['status']=='FRESH_GLOBAL_VARIATIONAL_SADDLE_CONNECTION_PASS','fresh original connection not complete')
 summaries=[];allvalues={};constants={}
 for name,side,clean in [('baseline','baseline_exact.hex','prerequisites/baseline.json'),('saddle','saddle_exact.hex','prerequisites/saddle.json'),('barrier','barrier_exact.hex','prerequisites/barrier.json'),('join_points','join_exact.hex','connections/join_points.json'),('leg_return','leg_return.exact.hex','connections/leg_return.json'),('leg_annihilation','leg_annihilation.exact.hex','connections/leg_annihilation.json'),('full_annihilation','full_annihilation.exact.hex','connections/full_annihilation.json')]:
  report,values,raw=summary(name,out/side,CLEAN/clean);summaries.append(report);allvalues[name]=values;constants[name]=raw
 # Enforce coverage of actual scalar consumers, independently of source schema.
 required={'baseline':['rho','delta','mu','branch_margin','c'],'saddle':['delta','negative_direction_bound','transverse_coercivity','rho','candidate_branch_margin','candidate_relative_to_reference/0','candidate_relative_to_reference/1'],'barrier':['marked_lower','longitudinal_upper','lateral_lower'],'leg_return':['path_relative_reference_upper'],'leg_annihilation':['path_relative_reference_upper']}
 for name,keys in required.items():
  for key in keys:need(key in allvalues[name],'missing consumed scalar '+name+':'+key)
 for k in range(2):
  for key in ('distance_to_saddle_candidate/1','signed_projection_on_normalized_q/0','signed_projection_on_normalized_q/1','relative_energy_to_reference/1'):need('points/'+str(k)+'/'+key in allvalues['join_points'],'missing join consumer')
 # Exact closure of compiled binary constant dependencies.
 b=allvalues['baseline'];s=allvalues['saddle'];cb=constants['barrier'];cp=constants['baseline']
 need(cp['probe_c']==b['c'] and cp['probe_rho']==b['rho'],'scalar probe arguments differ from previous exact baseline values')
 need(b['rho']<cb['compiled_rho'],'baseline rho exceeds compiled barrier radius')
 need(b['delta']<cb['compiled_delta'],'baseline delta exceeds compiled barrier delta')
 need(b['mu']>cb['compiled_mustar'],'baseline mu below compiled barrier coercivity')
 need(b['branch_margin']==cp['barrier_branch_limit'],'baseline branch margin differs from compiled barrier branch limit')
 for key in ('compiled_rho','compiled_delta','compiled_mustar'):need(cb[key]==cp[key],'compiled literal mismatch between probe and actual barrier program')
 gaps={'barrier_rho_minus_baseline_rho':str(cb['compiled_rho']-b['rho']),'barrier_delta_minus_baseline_delta':str(cb['compiled_delta']-b['delta']),'baseline_mu_minus_barrier_mustar':str(b['mu']-cb['compiled_mustar']),'baseline_and_barrier_branch_literal_identical':True}
 # Optimality overlay compares bare decimals to exact rational STRING constants.
 for name,raw,threshold,dr in [('minimum_rho',b['rho'],F('.017422'),'U'),('minimum_delta',b['delta'],F('.001166548'),'U'),('saddle_rho',s['rho'],F('.012825'),'U'),('saddle_positive_margin',s['transverse_coercivity'],F('.0719549999'),'L')]:
  gap=threshold-raw if dr=='U' else raw-threshold;need(gap>0,'exact unpadded overlay threshold fails '+name);gaps['optimality_'+name+'_strict_gap']=str(gap)
 # Supplementary node/segment tables are approximate decimal records; verify PAD
 # enclosure here, without relabeling their bare decimal endpoints as outward.
 table_checks=[]
 for name in ('leg_return','leg_annihilation','full_annihilation'):
  nodes={d['node']:d for d in [json.loads(x,parse_float=F) for x in (out/(name+'.nodes.jsonl')).read_text().splitlines()]};segs={d['segment']:d for d in [json.loads(x,parse_float=F) for x in (out/(name+'.segments.jsonl')).read_text().splitlines()]}
  for ext in ('nodes.jsonl','segments.jsonl'):need((out/(name+'.'+ext)).read_bytes()==(CLEAN/'connections'/(name+'.'+ext)).read_bytes(),'diagnostic table changed '+name)
  for line in (out/(name+'.exact.hex')).read_text().splitlines():
   v=line.split()
   if v[0]=='NODE':
    k=int(v[1]);q=nodes.pop(k);need(hf(v[2])<=hf(v[3]),'unordered node energy interval')
    for j,dr in ((0,'L'),(1,'U')):table_checks.append(one_bound(v[2+j],q['relative_energy'][j],dr,f'{name}:node{k}:energy{j}'))
    table_checks.append(one_bound(v[4],q['gradient_norm_upper'],'U',f'{name}:node{k}:norm'))
   elif v[0]=='SEG':
    k=int(v[1]);q=segs.pop(k);table_checks.append(one_bound(v[2],q['curvature_error_upper'],'U',f'{name}:segment{k}:curvature'));table_checks.append(one_bound(v[3],q['relative_energy_upper'],'U',f'{name}:segment{k}:energy'))
  need(not nodes and not segs,'incomplete exact path table coverage')
 sourcefiles=[p for p in ROOT.rglob('*') if p.is_file() and p.suffix in ('.cpp','.py','.patch') and 'results' not in p.parts and '__pycache__' not in p.parts]
 result={'status':'CONNECTION_FIXED_EXPORTS_EXACT_COMPARISON_PASS','scope':'fixed scientific inputs, original arithmetic plus finite/order guards and exact-hex diagnostics; not a second interval implementation','all_seven_scientific_summary_json_byte_identical_to_fresh_original':True,'all_six_path_diagnostic_tables_byte_identical_to_fresh_original':True,'summary_scalar_checks':sum(len(x['checks']) for x in summaries),'path_table_scalar_checks':len(table_checks),'summaries':summaries,'path_table_checks':table_checks,'compiled_and_unpadded_dependencies':gaps,'baseline_guard_reuse':reuse,'selected_saddle_guard':selected,'fixed_input_bindings':inputbindings,'flags':json.loads((ROOT/'instrumentation.json').read_text())[0]['flags'],'frozen_scientific_manifest_sha256':sha(SRC/'MANIFEST.sha256'),'fresh_original_connection_certificate_sha256':sha(CLEAN/'certificate.json'),'audit_source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sourcefiles}}
 (out/'exact_export_report.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],result['summary_scalar_checks'],result['path_table_scalar_checks'])
if __name__=='__main__':main()
