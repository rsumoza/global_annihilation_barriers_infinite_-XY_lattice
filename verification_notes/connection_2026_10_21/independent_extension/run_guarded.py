#!/usr/bin/env python3
"""Run guarded extension modules only; never invoke the global supervisor."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import re
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parent
WORK=ROOT/'work'
OUT=ROOT/'results'
CHANGED=(42,43,44,53,55)

def need(ok,reason):
    if not ok:raise RuntimeError(reason)

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def run(command,cwd,output,log):
    with output.open('w') as out,log.open('w') as err:
        p=subprocess.run([str(x) for x in command],cwd=cwd,stdout=out,stderr=err,timeout=1200)
    need(p.returncode==0,'guarded module rejected: '+str(command)+'; '+str(log))
    counters=re.search(r'GUARD_DIAGNOSTIC finite_interval_constructions=(\d+) finite_scalar_calls=(\d+)',log.read_text())
    need(counters is not None,'missing completed guard counters')
    return {'interval_constructions':int(counters[1]),'dn_up_calls':int(counters[2])}

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=True);records=[]
    qold=WORK/'common/Q_old.raw';qnew=WORK/'common/q_negative.raw';hs=WORK/'common/saddle_h.raw';hb=WORK/'common/base_h.raw'
    for ref,name in ((0,'plus'),(1,'minus')):
        for stem,args in [
            ('projection',[qold,qnew,hs,hb,ref]),
            ('center_distances',[qold,hs,hb,WORK/'centers',ref]),
            ('adjoint',[1024,qold,qnew,hs,ref]),
        ]:
            label=stem+'_'+name;print('Guarded '+label,flush=True)
            t=time.monotonic();counters=run([WORK/('guarded_'+stem),*args],WORK,OUT/(label+'.json'),OUT/(label+'.log'))
            exact=OUT/(label+'.exact.hex');shutil.copyfile(WORK/'exact_exports.hex',exact)
            need(exact.stat().st_size>0,'empty exact module sidecar')
            records.append({'label':label,'kind':stem,'reflection':ref,'duration_seconds':time.monotonic()-t,'counters':counters,'json_sha256':sha(OUT/(label+'.json')),'exact_sidecar_sha256':sha(exact)})
    chartout=OUT/'enlarged_charts';chartout.mkdir(exist_ok=True)
    def chart(cid):
        label=f'c{cid:03}';print('Guarded chart '+label,flush=True)
        t=time.monotonic();cwd=WORK/'chart_work'/label;dest=chartout/label;dest.mkdir(exist_ok=True)
        counters=run([WORK/'guarded_chart',1],cwd,dest/'chart.json',dest/'chart.log')
        status=json.loads((dest/'chart.json').read_text());need(status['interior_pass'] and status['mu']>0,'chart prerequisites')
        for name in ('coefficients.txt','coefficients_exact.hex'):shutil.copyfile(cwd/'results'/name,dest/name)
        return {'label':label,'kind':'chart','duration_seconds':time.monotonic()-t,'counters':counters,'json_sha256':sha(dest/'chart.json'),'exact_sidecar_sha256':sha(dest/'coefficients_exact.hex')}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:records.extend(pool.map(chart,CHANGED))
    result={'status':'GUARDED_EXTENSION_FIXED_INPUTS_PASS','scope':'Original arithmetic formulas plus finite/order instrumentation; not an independent interval backend or global certificate','duration_seconds':time.monotonic()-start,'runs':records,'total_interval_constructions':sum(x['counters']['interval_constructions'] for x in records),'total_dn_up_calls':sum(x['counters']['dn_up_calls'] for x in records)}
    (OUT/'guarded_replay.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],flush=True)

if __name__=='__main__':main()
