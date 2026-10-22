#!/usr/bin/env python3
"""Pinned source audit and conditional scalar replay, not fresh C++ proofs."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent/'xy_global_connection'
D=ROOT/'distribution';P=ROOT/'results_provided'
MPIN='ba100e9b7b8763bcd9ce1bcb7ee1a50a1f76669d2be92966ecbe6f6208feb27c'
SPIN='7a2f72741a64846029d6ebd79dd338ce15a587f48fb7fe71e6081af336040d7c'
INNERPIN='6a8b3f625b57032b76ca30f58c37bdd71115ea4a4bc682d0d7ac7c76b4ceb5a5'
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def require(ok,why):
    if not ok:raise RuntimeError(why)
def load(p):return json.loads(p.read_text(),parse_float=F)
def run(cmd,label):
    with (HERE/(label+'.stdout.log')).open('w') as out,(HERE/(label+'.stderr.log')).open('w') as err:
        r=subprocess.run(cmd,stdout=out,stderr=err)
    require(r.returncode==0,label+' rejected')

require(sha(D/'MANIFEST.sha256')==MPIN,'connection manifest pin')
require(sha(ROOT/'verify_global.py')==SPIN,'connection supervisor pin')
require(sha(D/'stage2/MANIFEST.sha256')==INNERPIN,'inner manifest pin')
records=[]
for line in (D/'MANIFEST.sha256').read_text().splitlines():
    digest,name=line.split('  ',1)
    require(name not in records,'duplicate member')
    require(sha(D/name)==digest,'member changed: '+name)
    records.append(name)
pins={'stage2/common/base_h.raw':'2a969db50023b0697c5f228d4b10325df7ea899bcdfe965fef6d02616f82adf3',
      'stage2/saddle/inputs/h.raw':'473f3ad414bdb6e59597b3b278628758d18e1563036075176df338d47c617cbb',
      'stage2/saddle/inputs/q.raw':'18152163150666b8108a35bf62da6f3b875875579fb4982a187d371abbcf45bd'}
for name,pin in pins.items():require(sha(D/name)==pin,'global reference binding')
nodebytes=256*256*8
legbind={}
for name,expected in [('leg_return',8),('leg_annihilation',19)]:
    cfg=list(map(int,(D/'connections'/(name+'.cfg')).read_text().split()))
    require(cfg==[expected,256,12],'leg config')
    p=D/'connections'/(name+'.z.raw')
    require(p.stat().st_size==expected*nodebytes,'leg raw length')
    with p.open('rb') as f:
        first=f.read(nodebytes);f.seek(-nodebytes,2);last=f.read(nodebytes)
    if name=='leg_return':require(first==bytes(nodebytes),'return reference is not exactly zero')
    join=last if name=='leg_return' else first
    legbind[name]={'nodes':expected,'extracted_join_node_sha256':hashlib.sha256(join).hexdigest(),
                   'join_node':'last' if name=='leg_return' else 'first'}

pre=P/'prerequisites';con=P/'connections'
rat=HERE/'rational_replayed.json';glob=HERE/'global_rational_replayed.json'
run([sys.executable,str(D/'stage2/math/postprocess_certificate.py'),str(pre/'saddle.json'),
     '--minimum-certificate',str(pre/'baseline.json'),'--output',str(rat)],'local_replay')
require(load(rat)==load(pre/'rational.json'),'local rational output differs')
run([sys.executable,str(D/'postprocess_global.py'),str(pre),str(con),'--output',str(glob)],'global_replay')
require(load(glob)==load(P/'global_rational.json'),'global rational output differs')

negative=[]
with tempfile.TemporaryDirectory() as td:
    tp=Path(td)/'prerequisites';tq=Path(td)/'connections';tp.mkdir();tq.mkdir()
    originals={}
    for folder,source,names in [(tp,pre,['baseline.json','saddle.json','rational.json','barrier.json']),
                                (tq,con,['leg_return.json','leg_annihilation.json','join_points.json'])]:
        for name in names:originals[folder/name]=json.loads((source/name).read_text())
    controls=[('not_index_one',tp/'saddle.json',lambda x:x.update(morse_index=2)),
              ('missing_time_coverage',tq/'leg_return.json',lambda x:x.update(all_time_intervals_covered=False)),
              ('nonzero_lift',tq/'leg_annihilation.json',lambda x:x.update(compact_integer_lift_max=1)),
              ('charged_endpoint',tq/'leg_annihilation.json',lambda x:x.update(endpoint_charge_count=1)),
              ('leg_above_saddle',tq/'leg_return.json',lambda x:x.update(path_relative_reference_upper=.006)),
              ('bad_join_direction',tq/'join_points.json',lambda x:x['points'][0].update(signed_projection_on_normalized_q=[0,0])),
              ('bad_Hermite_energy',tq/'join_points.json',lambda x:x['points'][0].update(relative_energy_to_reference=[.01,.01]))]
    for label,path,mutate in controls:
        for p,obj in originals.items():p.write_text(json.dumps(obj))
        obj=json.loads(path.read_text());mutate(obj);path.write_text(json.dumps(obj))
        out=Path(td)/'result.json'
        r=subprocess.run([sys.executable,str(D/'postprocess_global.py'),str(tp),str(tq),'--output',str(out)],capture_output=True,text=True)
        require(r.returncode!=0 and not out.exists(),'invalid connection accepted: '+label)
        negative.append({'case':label,'returncode':r.returncode,'message':r.stderr.strip()})

base=load(pre/'baseline.json');saddle=load(pre/'saddle.json');r=load(glob)
pad=F('1e-18')
consumed_saddle=[F(saddle[k]) for k in ['delta','negative_direction_bound','transverse_coercivity','rho','candidate_branch_margin']]
require(all(0<=x<1 for x in consumed_saddle),'saddle serialization magnitude')
join=load(con/'join_points.json')
joinvals=[F(v) for point in join['points'] for key in ['distance_to_saddle_candidate','signed_projection_on_normalized_q','relative_energy_to_reference'] for v in point[key]]
require(all(abs(x)<1 for x in joinvals),'join serialization magnitude')
thresholds={'minimum_rho':F('.017422')-base['rho']-pad,
            'minimum_delta':F('.001166548')-base['delta']-pad,
            'cap_rho':F('.012825')-F(saddle['rho'])-pad,
            'cap_candidate_coercivity':F(saddle['transverse_coercivity'])-pad-F('.0719549999')}
require(all(x>0 for x in thresholds.values()),'padded optimality contract fails')
levelslack=F('.00568947710965782931')-F(r['exact_saddle_relative_reference'][1]['exact'])
require(levelslack>0,'exact outer-level comparison')
summary={'status':'INDEPENDENT_CONNECTION_BINDINGS_AND_CONDITIONAL_SCALARS_PASS',
 'scope':'Source/raw bindings and replay of scalar consequences on archived local enclosures; not a fresh C++ CAP.',
 'manifest_sha256':MPIN,'supervisor_sha256':SPIN,'inner_manifest_sha256':INNERPIN,
 'scientific_members':len(records),'fixed_raw_bindings':pins,'actual_leg_endpoint_bindings':legbind,
 'local_rational_replay_matches':True,'global_rational_replay_matches':True,
 'padded_threshold_slacks':{k:str(x) for k,x in thresholds.items()},
 'outer_level_minus_padded_exact_saddle_upper':str(levelslack),
 'Hermite_margins':{v['label']:v['Hermite_strict_margin']['exact'] for v in r['join_points']},
 'face_attachment_margins':{v['label']:v['face_attachment_margin']['exact'] for v in r['join_points']},
 'leg_margins':{k:v['strict_margin_below_saddle']['exact'] for k,v in r['global_legs'].items()},
 'negative_controls':negative}
(HERE/'audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
