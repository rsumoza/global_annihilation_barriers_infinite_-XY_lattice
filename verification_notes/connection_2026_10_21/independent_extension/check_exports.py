#!/usr/bin/env python3
"""Exact fixed-export audit and scientific byte comparison for the extension.

This checks serialization of bounds from a guarded original engine. It does
not independently calculate the original interval bounds or analytic tails.
"""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math

ROOT=Path(__file__).resolve().parent
CHANGED=(42,43,44,53,55)
PAD=F('1e-18')
LEVEL=F('.00568947710965782931')

def need(ok,reason):
    if not ok:raise RuntimeError(reason)

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def load(path):return json.loads(path.read_text(),parse_float=F)

def hf(value):
    sign=-1 if value.startswith('-') else 1
    value=value.lstrip('+-');mant,ex=value.lower().split('p');need(mant.startswith('0x'),'invalid hex dyadic')
    a,_,b=mant[2:].partition('.');q=F(sign*int(a+b,16),16**len(b));ex=int(ex)
    return q*(2**ex if ex>=0 else F(1,2**(-ex)))

def interval_records(obj,path=()):
    if isinstance(obj,list):
        if len(obj)==2 and all(isinstance(x,str) for x in obj):
            try:
                yield path,tuple(map(F,obj))
                return
            except ValueError:pass
        for i,v in enumerate(obj):yield from interval_records(v,path+(str(i),))
    elif isinstance(obj,dict):
        for key,v in obj.items():yield from interval_records(v,path+(key,))

def compare_interval(raw,printed,label):
    a,b=map(hf,raw);p,q=map(F,printed)
    need(a<=b,'binary interval order: '+label)
    need(p<=a and q>=b,'decimal export direction: '+label)

def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--clean',type=Path,default=ROOT.parent/'full_fresh')
    parser.add_argument('--separator-results',type=Path)
    parser.add_argument('--connection-results',type=Path)
    parser.add_argument('--label',choices=('fresh','archived'),default='fresh')
    parser.add_argument('--output',type=Path,default=ROOT/'results/exact_exports_report.json')
    args=parser.parse_args();clean=args.clean.resolve();out=ROOT/'results'
    sep=args.separator_results or clean/'fresh_separator';con=args.connection_results or clean/'fresh_connection'
    need(load(out/'guarded_replay.json')['status']=='GUARDED_EXTENSION_FIXED_INPUTS_PASS','guarded replay incomplete')
    modules=[];all_constants={};bindings=[];interval_count=0
    for ref,name in ((0,'plus'),(1,'minus')):
        for stem in ('projection','center_distances','adjoint'):
            label=stem+'_'+name;jsonpath=out/(label+'.json');exact=out/(label+'.exact.hex')
            need(jsonpath.read_bytes()==(clean/(label+'.json')).read_bytes(),'scientific JSON changed: '+label)
            parsed=load(jsonpath);printed=list(interval_records(parsed));raw=[];constants={}
            for line in exact.read_text().splitlines():
                v=line.split()
                if v[0]=='CONST':need(v[1] not in constants,'duplicate constant');constants[v[1]]=hf(v[2])
                elif v[0]=='I':need(int(v[1])==len(raw),'exact interval record ledger');raw.append(v[2:4])
                else:raise RuntimeError('unknown exact export tag')
            need(len(raw)==len(printed),'interval export record count differs: '+label)
            for i,(binary,(path,decimal)) in enumerate(zip(raw,printed)):
                compare_interval(binary,decimal,label+':'+'.'.join(path))
            all_constants[label]=constants;interval_count+=len(raw)
            modules.append({'label':label,'interval_records_checked':len(raw),'scientific_bytes_identical':True})
            bindings.append({'label':label,'scientific_json_sha256':sha(jsonpath),'exact_sidecar_sha256':sha(exact)})
    chart_interval_count=0;chart_scalars=0;chart_headers=[]
    for cid in CHANGED:
        stem=f'c{cid:03}';d=out/'enlarged_charts'/stem;old=clean/'enlarged_charts'/stem
        for filename in ('chart.json','coefficients.txt'):
            need((d/filename).read_bytes()==(old/filename).read_bytes(),'chart scientific bytes changed: '+stem+'/'+filename)
        rows={};header={}
        for line in (d/'coefficients.txt').read_text().splitlines():
            v=line.split();tag=v[0]
            if tag=='HEADER':header=dict(zip(('mu','delta','remE','remF'),v[4:8]))
            elif tag=='DOMAIN_HEX':continue
            elif tag in ('E','B'):rows[tuple(v[:3])]=v[3:5]
            elif tag=='A':rows[tuple(v[:2])]=v[2:4]
            elif tag=='Y2':rows[tuple(v[:1])]=v[1:3]
            else:raise RuntimeError('unknown coefficient record')
        binary_header={}
        for line in (d/'coefficients_exact.hex').read_text().splitlines():
            v=line.split();tag=v[0]
            if tag=='HEADER':
                key=v[1];raw=hf(v[2]);printed=F(header[key]);need(key not in binary_header,'duplicate exact header')
                need(printed<=raw if key=='mu' else printed>=raw,'header direction differs: '+stem+'/'+key)
                binary_header[key]=raw;chart_scalars+=1
            else:
                key=tuple(v[:-2]);need(key in rows,'exact/export coefficient ledger');compare_interval(v[-2:],rows.pop(key),stem+'/'+str(key));chart_interval_count+=1
        need(not rows and set(binary_header)=={'mu','delta','remE','remF'},'incomplete exact chart record')
        chart_headers.append(binary_header)
        bindings.append({'label':stem,'scientific_json_sha256':sha(d/'chart.json'),'scientific_coefficients_sha256':sha(d/'coefficients.txt'),'exact_sidecar_sha256':sha(d/'coefficients_exact.hex')})
    # Literal dependencies are checked using their ACTUAL compiled dyadic values.
    baseline=load(sep/'baseline.json');saddle=load(con/'prerequisites/saddle.json')
    rho=all_constants['projection_plus']['rho_safe'];m=all_constants['projection_plus']['m_safe'];candidate=all_constants['projection_plus']['candidate_complement']
    need(all_constants['projection_minus']==all_constants['projection_plus'],'reflection compiled constants differ')
    need(all_constants['adjoint_plus']['rho_safe']==rho==all_constants['adjoint_minus']['rho_safe'],'adjoint/projection radius differs')
    need(rho>saddle['rho']+PAD,'compiled safe rho does not exceed saddle radius')
    need(candidate<saddle['transverse_coercivity']-PAD,'compiled complement exceeds certified candidate lower bound')
    scale=1<<256;sqrt2_lower=F(math.isqrt(2*scale*scale),scale)
    m_slack=candidate-rho/sqrt2_lower-m
    need(m_slack>0,'compiled m is not safely below compressed complement')
    deltas={h['delta'] for h in chart_headers};need(len(deltas)==1,'chart residual literals differ')
    delta=deltas.pop();delta_slack=delta-baseline['delta']-PAD;need(delta_slack>0,'compiled chart residual dependency fails')
    rational=load(con/'global_rational.json')
    exact_saddle_upper=F(rational['exact_saddle_relative_reference'][1]['exact'])
    level_slack=LEVEL-exact_saddle_upper;need(level_slack>0,'outer LEVEL is not above exact saddle upper bound')
    source=load(ROOT/'instrumentation.json')
    for filename,record in source['source_hashes'].items():need(sha(ROOT/'instrumented'/filename)==record['guarded'],'guarded audit source changed')
    result={
        'status':'GUARDED_EXTENSION_EXACT_EXPORTS_'+args.label.upper()+'_COMPARISON_PASS',
        'scope':'Fixed extension execution, original arithmetic with finite/order guards; interval exports checked against exact binary values; not an independent interval backend',
        'comparison_target':args.label,'modules':modules,'charts':len(CHANGED),'module_interval_records':interval_count,'chart_interval_records':chart_interval_count,'chart_one_sided_scalars':chart_scalars,
        'all_scientific_outputs_byte_identical':True,'literal_dependencies':{'compiled_rho':str(rho),'rho_minus_padded_saddle_radius':str(rho-saddle['rho']-PAD),'compiled_candidate_complement':str(candidate),'padded_candidate_minus_compiled':str(saddle['transverse_coercivity']-PAD-candidate),'compiled_m':str(m),'compressed_lower_minus_compiled_m':str(m_slack),'compiled_chart_delta':str(delta),'compiled_delta_minus_padded_baseline_delta':str(delta_slack)},
        'tight_energy_threshold':{'LEVEL':str(LEVEL),'saddle_upper_exact':str(exact_saddle_upper),'positive_exact_slack':str(level_slack),'display_as_float_for_orientation_only':float(level_slack)},
        'compiler_flags':source['flags'],'guarded_source_sha256':{k:v['guarded'] for k,v in source['source_hashes'].items()},'bindings':bindings,
        'input_bindings_sha256':sha(ROOT/'input_bindings.json'),'guarded_replay_sha256':sha(out/'guarded_replay.json'),
    }
    for filename,field in [(clean/'certificate.json','clean_full_certificate_sha256'),(sep/'certificate.json','separator_certificate_sha256'),(con/'certificate.json','connection_certificate_sha256')]:
        if filename.is_file():result[field]=sha(filename)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('bindings','literal_dependencies','guarded_source_sha256')},indent=2))

if __name__=='__main__':main()
