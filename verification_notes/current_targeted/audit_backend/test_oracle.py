from pathlib import Path
from fractions import Fraction as F
import subprocess,random,json,math
from decimal import Decimal as D, getcontext
getcontext().prec=160

def dec_atan(x):
 t=x/(1+(1+x*x).sqrt());v=t;total=D(0)
 for k in range(500):
  total+=v/(2*k+1);v=-v*t*t
  if abs(v)<D(10)**-150:break
 return 2*total

def dec_sin(x):
 v=x;total=v
 for k in range(1,500):
  v=-v*x*x/((2*k)*(2*k+1));total+=v
  if abs(v)<D(10)**-150:break
 return total

def dec_cos(x):
 v=D(1);total=v
 for k in range(1,500):
  v=-v*x*x/((2*k-1)*(2*k));total+=v
  if abs(v)<D(10)**-150:break
 return total
PI=16*dec_atan(D(1)/5)-4*dec_atan(D(1)/239)
root=Path(__file__).resolve().parent
exe=root/'engine_oracle'
def mh(s):
 s=s.strip().lower()
 sign=-1 if s.startswith('-') else 1
 s=s.lstrip('+-')
 if s in ('inf','nan'):return D(s)*sign
 mant,e=s.split('p');e=int(e)
 a,b=mant[2:].split('.') if '.' in mant else (mant[2:],'')
 return sign*D(int(a+b,16))*D(2)**(e-4*len(b))
def mf(x):return mh(x.hex())
rng=random.Random(20261020)
cases=[]
for op in ('add','sub','mul','div'):
 for _ in range(120):
  x=sorted([rng.uniform(-10,10),rng.uniform(-10,10)])
  y=sorted([rng.uniform(.01,10),rng.uniform(.01,10)])
  cases.append((op,*x,*y))
for op in ('sin','cos'):
 for _ in range(180):
  a=rng.uniform(-7.999,7.999);b=a if rng.randrange(2) else min(7.999,a+rng.uniform(0,1))
  cases.append((op,a,b,0.,0.))
 cases.extend([(op,-7.99,7.99,0.,0.),(op,0.,0.,0.,0.),(op,.25,.25,0.,0.)])
for _ in range(160):
 a=rng.uniform(-1.001,1.001);b=min(1.001,a+rng.uniform(0,.0001))
 cases.append(('atan',a,b,0.,0.))
for _ in range(160):
 a=2.**rng.randint(-900,900)*rng.uniform(1,2);b=a*rng.uniform(1,2)
 cases.append(('sqrt',a,b,0.,0.))
text=''.join(op+' '+' '.join(x.hex() for x in (a,b,c,d))+'\n' for op,a,b,c,d in cases)
p=subprocess.run([str(exe),'--oracle'],input=text,text=True,capture_output=True)
assert p.returncode==0,p.stderr
lines=p.stdout.splitlines();assert len(lines)==len(cases)
fail=[]
for i,((op,a,b,c,d),line) in enumerate(zip(cases,lines)):
 lo,hi=map(mh,line.split());a,b,c,d=map(mf,(a,b,c,d))
 if op=='add':v=[a+c,b+d]
 elif op=='sub':v=[a-d,b-c]
 elif op=='mul':v=[a*c,a*d,b*c,b*d]
 elif op=='div':v=[a/c,a/d,b/c,b/d]
 elif op=='sqrt':v=[a.sqrt(),b.sqrt()]
 elif op=='atan':v=[dec_atan(a),dec_atan(b)]
 else:
  fn=dec_sin if op=='sin' else dec_cos;v=[fn(a),fn(b)]
  shift=PI/2 if op=='sin' else 0
  for k in range(-10,11):
   t=shift+k*PI
   if a<=t<=b:v.append(fn(t))
 if not(lo<=min(v) and max(v)<=hi):fail.append({'i':i,'op':op,'got':line,'reference':[str(min(v)),str(max(v))]})
# Exact rational checks on published fixed Taylor remainders.
tails={
 'cos_small_degree_16':F(1,4)**18/math.factorial(18)<F(1,2**80),
 'sin_small_degree_17':F(1,4)**19/math.factorial(19)<F(1,2**80),
 'atan_36_terms':F(43,100)**73/73<F(1,2**90),
 'arcsin_principal_remainder_geometric':F(1,10)**11/(1-F(1,10)**2)<2*F(1,10)**11,
}
rejects={}
for name,row in [('trig_outside_domain','sin 0x1p+3 0x1p+3 0x0p+0 0x0p+0\n'),('atan_reduction_wide','atan -0x1p+0 0x1p+0 0x0p+0 0x0p+0\n'),('root_negative','sqrt -0x1p+0 -0x1p+0 0x0p+0 0x0p+0\n'),('nan_input','add nan nan 0x0p+0 0x0p+0\n')]:
 r=subprocess.run([str(exe),'--oracle'],input=row,text=True,capture_output=True)
 rejects[name]={'code':r.returncode,'stderr':r.stderr.strip()}
report={'scope':'Finite independent high-precision diagnostic tests; not a proof of the entire interval contract. Inputs are exact binary64 dyadics converted exactly to long double. Decimal comparison at 160 digits is diagnostic, not a certified oracle. Taylor inequalities below use exact fractions.','cases':len(cases),'failures':fail,'exact_rational_tail_inequalities':tails,'domain_rejections':rejects,'independent_reference':'Python decimal, precision160; independently coded Taylor/half-angle formulas'}
(root/'oracle_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
assert not fail and all(tails.values())
