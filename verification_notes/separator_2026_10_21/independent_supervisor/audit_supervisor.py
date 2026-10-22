#!/usr/bin/env python3
"""Independent exact ledger checks and conditional downstream replay.

Uses archived local enclosures, not a replacement for a fresh C++ run.
The source package is read only. All new outputs go beside this script.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, subprocess, sys, tempfile, shutil

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent/'xy_global_separator'
D=ROOT/'distribution'
P=ROOT/'results_provided'
def require(x, why):
    if not x: raise RuntimeError(why)
def sha(p): return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def load(p): return json.loads(p.read_text(),parse_float=F)

manifest_pin='991c744eab7e812d629662a0c0bc5cbc177e9dd2e9daa67caeefab7a814f941b'
supervisor_pin='d4646d988e538075b51d0c0bdd2d970976d1ccb53156ad52d414ec87f5c0095c'
require(sha(D/'MANIFEST.sha256')==manifest_pin,'frozen scientific manifest differs')
require(sha(ROOT/'verify_separator.py')==supervisor_pin,'frozen supervisor differs')
rows=[]
for line in (D/'MANIFEST.sha256').read_text().splitlines():
    digest,name=line.split('  ',1)
    require(sha(D/name)==digest,'manifest payload mismatch: '+name)
    rows.append(name)
require(len(rows)==len(set(rows)),'duplicate scientific manifest entry')
cases=load(D/'cases.json');grid=load(D/'grid.json')
cfg=(D/'geometry/cfg').read_text().splitlines()
require(list(map(int,cfg[0].split()))==[128,109],'geometry dimensions')
require(len(cfg)==110,'geometry ledger length')
require([x['id'] for x in cases]==list(range(109)),'case ids')
require([x['id'] for x in grid['centers']]==list(range(109)),'grid ids')
for case,gc,cr in zip(cases,grid['centers'],cfg[1:]):
    cid,gi,gj,name=cr.split()
    require(case['id']==gc['id']==int(cid),'id binding')
    require(case['grid_index']==gc['grid_index']==[int(gi),int(gj)],'grid binding')
    require(case['center']==gc['file']==name,'center bytes binding')
    vals=(D/'cases'/case['stem']/'cfg').read_text().split()
    require(int(vals[0])==128,'chart dimension')
    for text in vals[1:]:
        q=F.from_float(float.fromhex(text))
        require(q>0,'positive exact binary64 config')
    for text in vals[1:3]:
        require(float.fromhex(text).hex()==text,'canonical binary64 radius')
columns=load(D/'common/green_columns.json')
require(len(columns['edges'])==len({tuple(x) for x in columns['edges']}),'duplicate Green edge')
require((D/'common/green_columns.raw').stat().st_size==32768*len(columns['edges']),'Green data size')

# Every baseline field consumed by postprocess has an explicit rational pad.
baseline=load(P/'baseline.json')
consumed={'delta':baseline['delta'],'rho':baseline['rho'],'c':baseline['c'],
          'energy_candidate.lower':baseline['energy_candidate'][0]}
require(all(abs(v)<100 for v in consumed.values()),'22-digit field magnitude')
decimal_error=F(5,10**20) # bound for nearest 22 significant digits with |x|<100
pad=F('1e-18')
require(decimal_error<pad,'pad does not dominate decimal roundoff')
require(baseline['delta']+pad<F('.001166548'),'residual threshold')
require(baseline['rho']+pad<F('.017422'),'minimum radius threshold')

cmd=[sys.executable,str(D/'check_atlas.py'),'--manifest',str(P/'atlas_manifest.json'),
     '--proposals',str(P/'grid.json'),'--geometry',str(P/'geometry.json'),
     '--boundary-mode','faces','--max-depth','7','--level','5500000001/1000000000000',
     '--output',str(HERE/'boundary_replayed.json')]
with (HERE/'boundary_replay.stdout.log').open('w') as out,(HERE/'boundary_replay.stderr.log').open('w') as err:
    run=subprocess.run(cmd,stdout=out,stderr=err)
require(run.returncode==0,'conditional boundary replay failed')
fresh=load(HERE/'boundary_replayed.json');old=load(P/'boundary.json')
for obj in (fresh,old):
    for key in ('elapsed_seconds','input_hashes'):obj.pop(key,None)
    for chart in obj['charts']:chart.pop('path',None)
require(fresh==old,'conditional boundary differs from archived output')

negative=[]
with tempfile.TemporaryDirectory() as td:
    tmp=Path(td)
    for name in ('geometry.json','boundary.json','atlas_manifest.json'):
        (tmp/name).symlink_to((P/name).resolve())
    (tmp/'charts').symlink_to((P/'charts').resolve(),target_is_directory=True)
    bads=[('nan_delta','delta',float('nan')),('too_large_delta','delta',.001166548),
          ('too_large_rho','rho',.017422),('nonpositive_c','c',0)]
    original=json.loads((P/'baseline.json').read_text())
    for label,key,value in bads:
        b=dict(original);b[key]=value
        (tmp/'baseline.json').write_text(json.dumps(b))
        r=subprocess.run([sys.executable,str(D/'postprocess_separator.py'),str(tmp),
                          '--output',str(tmp/'result.json')],capture_output=True,text=True)
        require(r.returncode!=0,'bad baseline accepted: '+label)
        require(not (tmp/'result.json').exists(),'bad baseline wrote pass output: '+label)
        negative.append({'case':label,'returncode':r.returncode,'message':r.stderr.strip()})

result={'status':'INDEPENDENT_SUPERVISOR_AND_CONDITIONAL_REPLAY_PASS',
        'scope':'Pinned-input identities and conditional replay only; no fresh C++ local proof here.',
        'scientific_manifest_sha256':manifest_pin,'supervisor_sha256':supervisor_pin,
        'scientific_manifest_members':len(rows),'bound_case_center_records':len(cases),
        'baseline_consumed_fields':{k:str(v) for k,v in consumed.items()},
        'baseline_decimal_roundoff_bound':str(decimal_error),'postprocess_pad':str(pad),
        'baseline_residual_threshold_slack':str(F('.001166548')-baseline['delta']-pad),
        'conditional_charts':len(fresh['charts']),
        'coordinate_face_evaluations':sum(x['coordinate_faces']['counts']['evaluated'] for x in fresh['charts']),
        'unresolved_radial_circle_cells':sum(x['counts']['unresolved'] for x in fresh['charts']),
        'unresolved_coordinate_face_cells':sum(x['coordinate_faces']['counts']['unresolved'] for x in fresh['charts']),
        'negative_controls':negative}
(HERE/'audit_summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
