from pathlib import Path
import json,shutil,hashlib,difflib,subprocess
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'xy_global_separator'/'distribution'
flags=['-std=c++17','-O2','-DNDEBUG','-fno-fast-math','-ffp-contract=off','-frounding-math']
records=[]
for name in ('baseline','chart','geometry'):
 d=ROOT/name;d.mkdir(exist_ok=True)
 s=(SRC/name/'engine.cpp').read_text();old=s
 s=s.replace('static R dn(R x){return std::nextafter(x,-std::numeric_limits<R>::infinity());}', '''static unsigned long long guarded_scalar_calls=0,guarded_interval_constructions=0;
static R dn(R x){guarded_scalar_calls++;need(std::isfinite(x),"nonfinite dn input");R y=std::nextafter(x,-std::numeric_limits<R>::infinity());need(std::isfinite(y),"nonfinite dn output");return y;}''')
 s=s.replace('static R up(R x){return std::nextafter(x,std::numeric_limits<R>::infinity());}', 'static R up(R x){guarded_scalar_calls++;need(std::isfinite(x),"nonfinite up input");R y=std::nextafter(x,std::numeric_limits<R>::infinity());need(std::isfinite(y),"nonfinite up output");return y;}')
 s=s.replace('struct I{R l,h;I(R a=0):l(a),h(a){} I(R a,R b):l(a),h(b){}};', '''struct I{R l,h;I(R a=0):l(a),h(a){guarded_interval_constructions++;need(std::isfinite(a),"nonfinite scalar interval");} I(R a,R b):l(a),h(b){guarded_interval_constructions++;need(std::isfinite(a)&&std::isfinite(b)&&a<=b,"invalid interval endpoints");}};
static void report_guards(){std::cerr<<"GUARD_DIAGNOSTIC finite_interval_constructions="<<guarded_interval_constructions<<" finite_scalar_calls="<<guarded_scalar_calls<<"\\n";}''')
 s=s.replace('std::vector<I> av(m*m),t(m*m);R invert=1e10L,pivot=1e10L;', '''need(A.size()==(size_t)m*m&&V.size()==(size_t)m*m,"congruence dimensions");
 for(const I& x:A)need(std::isfinite(x.l)&&std::isfinite(x.h)&&x.l<=x.h,"nonfinite congruence input A");
 for(double x:V)need(std::isfinite(x),"nonfinite congruence input V");
 std::vector<I> av(m*m),t(m*m);R invert=1e10L,pivot=1e10L;''')
 if name=='baseline':
  s=s.replace('out<<std::setprecision(22)', 'std::cerr<<std::hexfloat<<"EXACT_BASELINE delta="<<delta<<" rho="<<rho<<" c="<<c<<" E_lower="<<initialE.l<<"\\n"<<std::defaultfloat;\n out<<std::setprecision(22)')
  s=s.replace('run(argv[1],argc==3?argv[2]:"");return 0;', 'run(argv[1],argc==3?argv[2]:"");report_guards();return 0;')
 (d/'engine.cpp').write_text(s)
 (d/'engine.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile=f'original/{name}/engine.cpp',tofile=f'instrumented/{name}/engine.cpp')))
 if name!='baseline':
  c=(SRC/name/'certify.cpp').read_text()
  c=c.replace('return 0;', 'report_guards();return 0;')
  (d/'certify.cpp').write_text(c)
 source='engine.cpp' if name=='baseline' else 'certify.cpp'
 subprocess.run(['g++',*flags,source,'-o','verify_guarded'],cwd=d,check=True)
 records.append({'program':name,'source_original_sha256':hashlib.sha256(old.encode()).hexdigest(),'source_instrumented_sha256':hashlib.sha256(s.encode()).hexdigest(),'flags':flags})
(ROOT/'instrumentation.json').write_text(json.dumps(records,indent=2)+'\n')
# Construct each fixed case independently using frozen proposal layout.
cases=json.loads((SRC/'cases.json').read_text());cols=json.loads((SRC/'common/green_columns.json').read_text());raw=(SRC/'common/green_columns.raw').read_bytes();index={tuple(e):k for k,e in enumerate(cols['edges'])}
for c in cases:
 src=SRC/'cases'/c['stem'];case=ROOT/'cases'/c['stem'];inp=case/'inputs';inp.mkdir(parents=True,exist_ok=True);(case/'results').mkdir(exist_ok=True)
 for n in ('cfg','edges','basis.raw','poly_green.raw'):
  dst=inp/n
  if not dst.exists():dst.hardlink_to(src/n)
 for s,dst in [(SRC/'common/h.raw',inp/'h.raw'),(SRC/'common/q.raw',inp/'q.raw'),(SRC/'centers'/c['center'],inp/'center.raw')]:
  if not dst.exists():dst.hardlink_to(s)
 selected=b''.join(raw[index[tuple(e)]*32768:(index[tuple(e)]+1)*32768] for e in c['green_edges'])
 
 if hashlib.sha256(selected).hexdigest()!=c['green_sha256']:raise RuntimeError('reconstructed Green identity differs')
 (inp/'green.raw').write_bytes(selected)
for p in (SRC/'baseline').iterdir():
 if p.is_file() and p.name!='engine.cpp':
  dst=ROOT/'baseline'/p.name
  if not dst.exists():dst.hardlink_to(p)
p=ROOT/'baseline/min.h'
if not p.exists():p.hardlink_to(SRC/'common/h.raw')
(ROOT/'cases.json').write_text(json.dumps(cases))
print('instrumented copies built and fixed inputs reconstructed')
