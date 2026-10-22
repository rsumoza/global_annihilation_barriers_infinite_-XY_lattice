#!/usr/bin/env python3
"""Focused tests of the ORIGINAL cap-checker's nested predicates.

Extracting the actual AST preserves branch inequalities and avoids editing the
original. These tests do not validate archived cap/chart/exterior inputs.
"""
from pathlib import Path
import ast, importlib.util, json, sys
from fractions import Fraction

ROOT=Path(__file__).resolve().parents[1]/'xy_optimality_extension'
spec=importlib.util.spec_from_file_location('geometry_test_atlas',ROOT/'extension/check_atlas.py')
b=importlib.util.module_from_spec(spec);sys.modules[spec.name]=b;spec.loader.exec_module(b)
I=b.I; Z=b.ZERO; O=b.ONE; T=b.TWO
src=ROOT/'extension/check_capped_separator.py'
tree=ast.parse(src.read_text())
main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
selected=[n for n in main.body if isinstance(n,ast.FunctionDef) and n.name in ('linear','outside_cut','root_cap_test')]
ns={'I':I,'G':[[O,Z],[Z,O]],'caps':[]}
exec(compile(ast.Module(body=selected,type_ignores=[]),str(src),'exec'),ns)
F=I.rational
base={'normal':[O,Z],'offset':Z,'candidate_a':[Z,Z], 'Ydist2':{0:Z},'rho':F('1/100'),'Rcap':F('33/200')}
second=dict(base,normal=[Z,O],offset=F(100))
ns['caps'][:]=[base,second]
checks={}
def record(name,condition):
 checks[name]=bool(condition)
 if not condition: raise AssertionError(name)

record('outside_cut_exact_tangency_retained',not ns['outside_cut'](Z,Z))
record('outside_cut_touching_positive_interval_retained',not ns['outside_cut'](I(0,F('1/100').hi),Z))
record('outside_cut_touching_negative_interval_retained',not ns['outside_cut'](I(F('-1/100').lo,0),Z))
record('outside_cut_straddling_interval_retained',not ns['outside_cut'](I(F('-1/100').lo,F('1/100').hi),Z))
record('outside_cut_strict_positive_excluded',ns['outside_cut'](F('1/100'),Z))

# Use the ORIGINAL plane-loop exclusion expression, not a restatement.
plane_if=next(n for n in ast.walk(main) if isinstance(n,ast.If) and n.lineno==111)
plane_expr=compile(ast.Expression(plane_if.test),str(src),'eval')
def plane_skips(x,y):
 env=dict(ns,cap=base,gx=x,gy=y,radius2=F(10000),ell=ns['linear'](base,x,y)-base['offset'])
 return eval(plane_expr,env)
record('plane_exact_zero_retained',not plane_skips(Z,Z))
record('plane_positive_touch_retained',not plane_skips(I(0,F('1/100').hi),Z))
record('plane_negative_touch_retained',not plane_skips(I(F('-1/100').lo,0),Z))
record('plane_strict_positive_excluded',plane_skips(F('1/100'),Z))
record('plane_strict_negative_excluded',plane_skips(F('-1/100'),Z))

chart={'mu':O,'R':F('1/5')}
ok,bound=ns['root_cap_test'](base,chart,0,Z,Z,F('1/20'))
record('fullfiber_interior_root_and_cap_accepted',ok)
cap_equal=dict(base,Rcap=F('3/50'))  # rho + beta/mu = .01+.05=.06
ok,bound=ns['root_cap_test'](cap_equal,chart,0,Z,Z,F('1/20'))
record('cap_radius_equality_rejected',not ok)
ok,bound=ns['root_cap_test'](base,chart,0,Z,Z,F('1/5'))
record('fiber_interior_root_equality_rejected',not ok)

# The ORIGINAL Taylor-positivity check uses Lambda=1/sqrt(2).
mu=F('1/25');rad=F('33/200')
record('displayed_cap_radius_passes_sharp_contract',(F(3)*T.sqrt()*mu-rad).lo>0)
record('displayed_cap_radius_not_justified_by_coarse_Lambda_one',(F(3)*mu-rad).hi<0)
record('rational_tau_above_Lambda',2*Fraction(71,100)**2>1)
record('rational_cap_margin_exact',Fraction(1,50)-Fraction(71,100)*Fraction(33,200)/6==Fraction(19,40000))

# Interval TightPoly's natural/mean intersection must contain the same actual
# polynomial for any selected real coefficients. Exhaust selected exact dyadic
# coefficient and point values; this tests the implementation's dependency use.
rows=[{'i':0,'j':0,'lower':'-1/4','upper':'1/4'},
      {'i':1,'j':1,'lower':'-3/2','upper':'5/4'},
      {'i':2,'j':0,'lower':'-1','upper':'2'},
      {'i':0,'j':3,'lower':'-1/2','upper':'1/2'}]
p=b.TightPoly(rows); x=I(F('-3/4').lo,F('1/2').hi);y=I(F('-1/2').lo,F('3/4').hi)
v=p.eval(x,y)
for c0 in map(Fraction,['-1/4','0','1/4']):
 for c11 in map(Fraction,['-3/2','0','5/4']):
  for c20 in map(Fraction,['-1','0','2']):
   for c03 in map(Fraction,['-1/2','0','1/2']):
    for xx in map(Fraction,['-3/4','0','1/2']):
     for yy in map(Fraction,['-1/2','0','3/4']):
      exact=c0+c11*xx*yy+c20*xx**2+c03*yy**3
      assert Fraction(v.lo,b.UNIT)<=exact<=Fraction(v.hi,b.UNIT)
record('TightPoly_selected_fixed_coefficients_contained',True)

out={'status':'FOCUSED_GEOMETRY_CONTRACT_TESTS_PASS','checks':checks,'checked_source':str(src),
     'scope':'Original nested predicates and selected interval examples only; not a fresh CAP.',
     'source_sha256':__import__('hashlib').sha256(src.read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
