#!/usr/bin/env python3
"""Independent dyadic check of normal/adjoint coordinates for both reflections.

No supplied interval engine is imported. Products use the degree-four
Laplacian of zero extensions; a/sqrt(d) inclusions are reduced to exact
rational sign and square comparisons. Does not check cap coercivity/energies.
"""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import struct

ROOT=Path(__file__).resolve().parent
QOLD='7df80c3db341a9ece19c77939ad1c1ba857ac8dbb2b51d4a1435f26ee88bf7b0'
QNEW='18152163150666b8108a35bf62da6f3b875875579fb4982a187d371abbcf45bd'

def need(ok,reason):
    if not ok:raise RuntimeError(reason)

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def decode(path,count,pin):
    data=path.read_bytes();need(len(data)==8*count and hashlib.sha256(data).hexdigest()==pin,'raw identity differs')
    values=struct.unpack('<'+'d'*count,data);need(all(math.isfinite(v) for v in values),'nonfinite raw')
    return [v.as_integer_ratio() for v in values]

def lap(v):
    n=64;out=[]
    for i in range(n):
        for j in range(n):
            k=i*n+j;s=4*v[k]
            if i:s-=v[k-n]
            if i+1<n:s-=v[k+n]
            if j:s-=v[k-1]
            if j+1<n:s-=v[k+1]
            out.append(s)
    return out

def dot(a,b):return sum(x*y for x,y in zip(a,b))

def radical_enclosure(a,d,bounds,label):
    """Check lower<=a/sqrt(d)<=upper with no approximate square root."""
    lower,upper=map(F,bounds);need(d>0 and lower<=upper,'domain/order '+label)
    aa=a*a
    if a>=0:
        need(lower<=0 or lower*lower*d<=aa,'radical lower '+label)
        need(upper>=0 and upper*upper*d>=aa,'radical upper '+label)
    else:
        need(lower<=0 and lower*lower*d>=aa,'negative radical lower '+label)
        need(upper>=0 or upper*upper*d<=aa,'negative radical upper '+label)

def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--q-old',type=Path,default=ROOT/'work/common/Q_old.raw')
    parser.add_argument('--q-new',type=Path,default=ROOT/'work/common/q_negative.raw')
    parser.add_argument('--results',type=Path,default=ROOT/'results')
    parser.add_argument('--output',type=Path,default=ROOT/'results/exact_normal_report.json')
    args=parser.parse_args()
    old=decode(args.q_old,8192,QOLD);new=decode(args.q_new,4096,QNEW)
    denominator=max(d for n,d in old+new);need(denominator&(denominator-1)==0,'non-dyadic denominator')
    Q=[[n*(denominator//d) for n,d in old[k*4096:(k+1)*4096]] for k in range(2)]
    q=[n*(denominator//d) for n,d in new]
    for v in Q+[q]:
        need(all(v[j]==v[63*64+j]==v[j*64]==v[j*64+63]==0 for j in range(64)),'nonzero boundary')
    lQ=[lap(v) for v in Q];scale=denominator*denominator
    G=[[F(dot(Q[k],lQ[l]),scale) for l in range(2)] for k in range(2)]
    need(G[0][1]==G[1][0],'exact Gram asymmetry')
    determinant=G[0][0]*G[1][1]-G[0][1]*G[1][0];need(G[0][0]>0 and determinant>0,'Gram not positive')
    Gi=[[G[1][1]/determinant,-G[0][1]/determinant],[-G[1][0]/determinant,G[0][0]/determinant]]
    records=[];bindings=[]
    for ref,name in ((0,'plus'),(1,'minus')):
        qr=q if ref==0 else [q[(63-i)*64+j] for i in range(64) for j in range(64)]
        g=[F(dot(qr,v),scale) for v in lQ]
        u=[sum(Gi[k][l]*g[l] for l in range(2)) for k in range(2)]
        d=sum(g[k]*u[k] for k in range(2));need(d>0,'zero projected norm')
        # n=(-g)/sqrt(d), c=(-Gi*g)/sqrt(d). Verify G*c=n
        # by equality of the rational numerators, not interval overlap.
        need(all(sum(G[k][l]*(-u[l]) for l in range(2))==-g[k] for k in range(2)),'exact normal identity fails')
        need(sum(u[k]*G[k][l]*u[l] for k in range(2) for l in range(2))==d,'exact unit normal identity fails')
        projection_path=args.results/('projection_'+name+'.json');adjoint_path=args.results/('adjoint_'+name+'.json')
        projection=json.loads(projection_path.read_text());adjoint=json.loads(adjoint_path.read_text())
        need(projection['reflected_x']==adjoint['reflected_x']==ref,'reflection output differs')
        for k in range(2):
            radical_enclosure(-g[k],d,projection['normal_coordinates'][k],name+':normal'+str(k))
            radical_enclosure(-u[k],d,adjoint['r_coefficients'][k],name+':adjoint'+str(k))
        records.append({'reflection':ref,'normal_coordinate_enclosures':2,'adjoint_coordinate_enclosures':2,'exact_Gc_equals_normal_proved_by_rational_numerators':True,'exact_r_unit_energy_norm_proved':True,'g':list(map(str,g)),'projected_norm_squared':str(d)})
        bindings.append({'reflection':ref,'projection_sha256':sha(projection_path),'adjoint_sha256':sha(adjoint_path)})
    result={'status':'INDEPENDENT_EXACT_NORMAL_ADJOINT_COORDINATES_PASS','scope':'Common fixed Gram/normal/adjoint coordinate identity and 8 interval inclusions only; excludes energy, coercivity, Poisson residual and global topology','method':'Exact binary64 dyadics, degree-four zero-extension Laplacian, rational Gram inverse and sign-aware squared radical comparisons','raw_Q_old_sha256':QOLD,'raw_Q_new_sha256':QNEW,'checker_sha256':sha(Path(__file__)),'common_denominator_log2':denominator.bit_length()-1,'enclosures_checked':8,'records':records,'output_bindings':bindings}
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(result['status'])

if __name__=='__main__':main()
