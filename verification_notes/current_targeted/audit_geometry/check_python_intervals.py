"""Adversarial finite exact-arithmetic checks; not a formal interval proof."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util,json,random,sys
root=Path(__file__).resolve().parent
source=root.parent/'xy_optimality_extension/extension/check_atlas.py'
spec=importlib.util.spec_from_file_location('audited_interval_module',source)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
rng=random.Random(20261020)
checks=0

def enclosure(iv,lo,hi):
 global checks
 assert Q(iv.lo,m.UNIT)<=lo<=hi<=Q(iv.hi,m.UNIT)
 checks+=1

def sample():
 x=Q(rng.randint(-10000,10000),rng.randint(1,9999)); y=Q(rng.randint(-10000,10000),rng.randint(1,9999))
 lo,hi=sorted((x,y));return lo,hi,m.I.endpoints(str(lo),str(hi))

for _ in range(4000):
 a,b,A=sample();c,d,B=sample()
 enclosure(A,a,b);enclosure(B,c,d)
 enclosure(A+B,a+c,b+d);enclosure(A-B,a-d,b-c)
 values=[a*c,a*d,b*c,b*d];enclosure(A*B,min(values),max(values))
 values=[a*a,b*b];enclosure(A.square(),0 if a<=0<=b else min(values),max(values))
 if B.lo>0 or B.hi<0:
  values=[a/c,a/d,b/c,b/d];enclosure(A/B,min(values),max(values))
 N=A.square();S=N.sqrt()
 assert S.lo>=0 and Q(S.lo,m.UNIT)**2<=Q(N.lo,m.UNIT)
 assert Q(S.hi,m.UNIT)**2>=Q(N.hi,m.UNIT)
 checks+=1

for _ in range(500):
 xlo,xhi,X=sample(); ylo,yhi,Y=sample()
 rows=[];exact=[]
 for i,j in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(2,2),(0,4)]:
  lo,hi,_=sample();rows.append({'i':i,'j':j,'lower':str(lo),'upper':str(hi)})
  exact.append((i,j,lo+(hi-lo)*Q(rng.randint(0,16),16)))
 p=m.TightPoly(rows);iv=p.eval(X,Y)
 for tx,ty in [(0,0),(1,1),(1,0),(0,1),(Q(1,2),Q(1,2)),(Q(1,3),Q(2,3))]:
  x=xlo+(xhi-xlo)*tx;y=ylo+(yhi-ylo)*ty
  val=sum(c*x**i*y**j for i,j,c in exact)
  enclosure(iv,val,val)

out={'scope':'Finite adversarial exact rational enclosure checks, not a formal proof or C++ audit','seed':20261020,'random_interval_trials':4000,'random_polynomial_trials':500,'total_checked_enclosures_or_sqrt_relations':checks,'status':'PASS'}
(root/'python_interval_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
