"""Fixed-input diagnostic: every exported bound is compared to its exact binary value.

No decimal endpoint is read via float. Hex long doubles are parsed as dyadic Fractions.
The implementation is diagnostic: it does not replace an independent interval engine.
"""
from pathlib import Path
from fractions import Fraction as F
import json,re,hashlib
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'xy_global_separator'/'distribution'
CLEAN=ROOT.parent/'fresh_results'
PAD=F('1e-18')
def need(ok,why):
 if not ok:raise RuntimeError(why)
def hf(s):
 sign=-1 if s.startswith('-') else 1;s=s.lstrip('+-');mant,e=s.lower().split('p');need(mant.startswith('0x'),'invalid hex float')
 mant=mant[2:];a,sep,b=mant.partition('.');digits=a+b
 q=F(sign*int(digits,16),16**len(b));ex=int(e)
 return q*(2**ex if ex>=0 else F(1,2**(-ex)))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(),parse_float=F)
def check_interval(raw,printed,label):
 lo,hi=map(hf,raw);p,q=map(F,printed)
 need(lo<=hi,label+': invalid binary interval');need(p<=lo and q>=hi,label+': export is not outward')
def main():
 chart_count=0;chart_intervals=0;chart_scalars=0;hex_hashes={};bytechecks=[]
 for c in json.loads((ROOT/'cases.json').read_text()):
  case=ROOT/'cases'/c['stem'];res=ROOT/'results'/c['stem'];exact=res/'coefficients_exact.hex';printed=res/'coefficients.txt'
  need(exact.is_file(),f'missing exact sidecar {exact}')
  rows={};header=None
  for line in printed.read_text().splitlines():
   v=line.split()
   if v[0]=='HEADER':header=dict(zip(('mu','delta','remE','remF'),v[4:8]))
   elif v[0]=='DOMAIN_HEX':continue
   elif v[0] in ('E','B'):rows[(v[0],v[1],v[2])]=v[3:5]
   elif v[0]=='A':rows[(v[0],v[1])]=v[2:4]
   elif v[0]=='Y2':rows[(v[0],)]=v[1:3]
   else:raise RuntimeError('unknown exported row')
  exact_count=0
  for line in exact.read_text().splitlines():
   v=line.split();tag=v[0]
   if tag=='HEADER':
    label=v[1];raw=hf(v[2]);p=F(header[label]);need(p<=raw if label=='mu' else p>=raw,c['stem']+' header direction '+label);chart_scalars+=1;continue
   if tag in ('E','B'):key=tuple(v[:3]);raw=v[3:5]
   elif tag=='A':key=tuple(v[:2]);raw=v[2:4]
   elif tag=='Y2':key=tuple(v[:1]);raw=v[1:3]
   else:raise RuntimeError('unknown exact row')
   check_interval(raw,rows.pop(key),c['stem']+' '+str(key));exact_count+=1;chart_intervals+=1
  need(not rows,c['stem']+' unmatched exported intervals')
  bytechecks.append({'id':c['id'],'coefficient_bytes_equal_clean':printed.read_bytes()==(CLEAN/'charts'/c['stem']/'coefficients.txt').read_bytes(),'stdout_bytes_equal_clean':(res/'chart.json').read_bytes()==(CLEAN/'charts'/c['stem']/'chart.json').read_bytes()})
  need(all(list(bytechecks[-1].values())[1:]),c['stem']+' scientific output changed')
  hex_hashes[c['stem']]=sha(exact);chart_count+=1
 # Geometry exports, including every neighbor endpoint and every center.
 exact=ROOT/'results/geometry_exact.hex'
 need(exact.is_file(),'missing exact geometry sidecar');g=load(ROOT/'results/geometry.json');rows={}
 for name in ('G','Gi'):
  for k in range(4):rows[(name,str(k))]=g[name][k]
 for c in g['centers']:
  for k in range(2):rows[('A',str(c['id']),str(k))]=c['a'][k]
  rows[('Y2',str(c['id']))]=c['Y_norm_squared']
 for d in g['neighbors']:rows[('D',str(d['i']),str(d['j']))]=d['distance']
 maxkeys={'neighbor':'maximum_neighbor_Y_distance_upper','Y':'maximum_center_Y_norm_upper','coord':'maximum_nominal_coordinate_error_upper'};geometry_intervals=0;geometry_scalars=0
 for line in exact.read_text().splitlines():
  v=line.split();tag=v[0]
  if tag=='MAX':need(F(g[maxkeys[v[1]]])>=hf(v[2]),'geometry maximum not outward');geometry_scalars+=1;continue
  key=tuple(v[:-2]);check_interval(v[-2:],rows.pop(key),'geometry '+str(key));geometry_intervals+=1
 need(not rows,'geometry missing exact rows');need((ROOT/'results/geometry.json').read_bytes()==(CLEAN/'geometry.json').read_bytes(),'geometry scientific output changed')
 # Four baseline scalars consumed by the exact separator postprocessor.
 text=(ROOT/'results/baseline.log').read_text();line=next(x for x in text.splitlines() if x.startswith('EXACT_BASELINE '));binary=dict(w.split('=',1) for w in line.split()[1:]);b=load(ROOT/'results/baseline.json');vals={'delta':b['delta'],'rho':b['rho'],'c':b['c'],'E_lower':b['energy_candidate'][0]};baseline=[]
 for k,printed in vals.items():
  raw=hf(binary[k]);need(abs(printed-raw)<PAD,'baseline PAD insufficient '+k)
  baseline.append({'name':k,'raw_hex':binary[k],'raw_exact':str(raw),'printed_exact':str(printed),'absolute_export_error_exact':str(abs(printed-raw)),'pad_minus_error_exact':str(PAD-abs(printed-raw))})
 need((ROOT/'results/baseline.json').read_bytes()==(CLEAN/'baseline.json').read_bytes(),'baseline scientific output changed')
 probe=(ROOT/'probe_original.txt').read_text();compiled_delta=hf(next(x.split()[1] for x in probe.splitlines() if x.startswith('CHART_DELTA_HEX ')));delta=hf(binary['delta']);need(delta<compiled_delta,'baseline does not satisfy chart compiled literal dependency')
 result={'status':'FIXED_EXPORTS_EXACT_COMPARISON_PASS','scope':'as-run diagnostics on fixed separator inputs; same arithmetic engine with guards','charts':chart_count,'chart_interval_records':chart_intervals,'chart_one_sided_scalars':chart_scalars,'geometry_interval_records':geometry_intervals,'geometry_one_sided_scalars':geometry_scalars,'all_scientific_outputs_byte_identical_to_clean':True,'baseline_consumed_scalars':baseline,'compiled_chart_delta_exact':str(compiled_delta),'compiled_chart_delta_minus_baseline_exact':str(compiled_delta-delta),'chart_sidecar_sha256':hex_hashes,'geometry_sidecar_sha256':sha(exact),'chart_bytechecks':bytechecks}
 result['compiler_flags']=json.loads((ROOT/'instrumentation.json').read_text())[0]['flags']
 result['frozen_scientific_manifest_sha256']=sha(SRC/'MANIFEST.sha256')
 result['clean_separator_certificate_sha256']=sha(CLEAN/'certificate.json')
 result['audit_source_sha256']={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'prepare_instrumented.py',ROOT/'add_exact_sidecars.py',ROOT/'run_instrumented.py',ROOT/'check_exact_exports.py',ROOT/'probe_original.cpp',ROOT/'probe_guarded.cpp',ROOT/'baseline/engine.cpp',ROOT/'chart/engine.cpp',ROOT/'chart/certify.cpp',ROOT/'geometry/engine.cpp',ROOT/'geometry/certify.cpp']}
 result['consumed_sidecar_bindings']=[{'id':c['id'],'exact_hex_sha256':sha(ROOT/'results'/c['stem']/'coefficients_exact.hex'),'coefficients_sha256':sha(ROOT/'results'/c['stem']/'coefficients.txt'),'chart_stdout_sha256':sha(ROOT/'results'/c['stem']/'chart.json')} for c in json.loads((ROOT/'cases.json').read_text())]
 result['geometry_stdout_sha256']=sha(ROOT/'results/geometry.json')
 result['baseline_stdout_sha256']=sha(ROOT/'results/baseline.json')
 result['baseline_exact_hex_sidecar_log_sha256']=sha(ROOT/'results/baseline.log')
 result['compiled_literal_probe_sha256']=sha(ROOT/'probe_original.txt')
 (ROOT/'results/exact_export_report.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],chart_count,chart_intervals,geometry_intervals)
if __name__=='__main__':main()
