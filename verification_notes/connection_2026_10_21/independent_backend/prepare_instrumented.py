from pathlib import Path
import subprocess,json,hashlib,difflib,os
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'xy_global_connection'/'distribution'
FLAGS=['-std=c++17','-O2','-DNDEBUG','-fno-fast-math','-ffp-contract=off','-frounding-math']
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def guard_engine(s):
 s=s.replace('static R dn(R x){return std::nextafter(x,-std::numeric_limits<R>::infinity());}', '''static unsigned long long guarded_scalar_calls=0,guarded_interval_constructions=0;
static R dn(R x){guarded_scalar_calls++;need(std::isfinite(x),"nonfinite dn input");R y=std::nextafter(x,-std::numeric_limits<R>::infinity());need(std::isfinite(y),"nonfinite dn output");return y;}''')
 s=s.replace('static R up(R x){return std::nextafter(x,std::numeric_limits<R>::infinity());}', 'static R up(R x){guarded_scalar_calls++;need(std::isfinite(x),"nonfinite up input");R y=std::nextafter(x,std::numeric_limits<R>::infinity());need(std::isfinite(y),"nonfinite up output");return y;}')
 s=s.replace('struct I{R l,h;I(R a=0):l(a),h(a){} I(R a,R b):l(a),h(b){}};', '''struct I{R l,h;I(R a=0):l(a),h(a){guarded_interval_constructions++;need(std::isfinite(a),"nonfinite scalar interval");} I(R a,R b):l(a),h(b){guarded_interval_constructions++;need(std::isfinite(a)&&std::isfinite(b)&&a<=b,"invalid interval endpoints");}};
static void report_guards(){std::cerr<<"GUARD_DIAGNOSTIC finite_interval_constructions="<<guarded_interval_constructions<<" finite_scalar_calls="<<guarded_scalar_calls<<"\\n";}''')
 s=s.replace('std::vector<I> av(m*m),t(m*m);R invert=1e10L,pivot=1e10L;', '''need(A.size()==(size_t)m*m&&V.size()==(size_t)m*m,"congruence dimensions");
 for(const I& x:A)need(std::isfinite(x.l)&&std::isfinite(x.h)&&x.l<=x.h,"nonfinite congruence A");
 for(double x:V)need(std::isfinite(x),"nonfinite congruence V");
 std::vector<I> av(m*m),t(m*m);R invert=1e10L,pivot=1e10L;''')
 return s

def sidecar(lines,name='results/exact.hex'):
 return 'std::ofstream ee("'+name+'");ee<<std::hexfloat;\n '+''.join('ee<<"'+label+' "<<'+expr+'<<"\\n";\n ' for label,expr in lines)+'ee.close();need((bool)ee,"exact sidecar output");\n '
records=[]
for label in ('stage2/baseline','stage2/saddle','stage2/barrier','connections'):
 src=SRC/label;dst=ROOT/label;dst.mkdir(parents=True,exist_ok=True);old=(src/'engine.cpp').read_text();s=guard_engine(old)
 if label.endswith('baseline'):
  dump=sidecar([('U delta','delta'),('U rho','rho'),('L mu','mu'),('L c','c'),('L branch_margin','margin'),('L energy_candidate/0','initialE.l'),('U energy_candidate/1','(initialE+tailE).h')])
  s=s.replace('out<<std::setprecision(22)',dump+'out<<std::setprecision(22)')
  s=s.replace('run(argv[1],argc==3?argv[2]:"");return 0;','run(argv[1],argc==3?argv[2]:"");report_guards();return 0;')
 (dst/'engine.cpp').write_text(s)
 (dst/'engine.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='original/'+label+'/engine.cpp',tofile='instrumented/'+label+'/engine.cpp')))
 if label.endswith('baseline'):programs=[('engine.cpp','verify_guarded')]
 elif label=='connections':programs=[('certify_path.cpp','path_guarded'),('certify_join_points.cpp','join_guarded')]
 else:programs=[('certify.cpp','verify_guarded')]
 for file,exe in programs:
  if file!='engine.cpp':
   old=(src/file).read_text();s=old
   if label.endswith('saddle'):
    dump=sidecar([('U delta','delta'),('L negative_direction_bound','negative'),('L transverse_coercivity','positive'),('U rho','rho'),('L candidate_branch_margin','margin'),('L candidate_relative_to_reference/0','relativeEnergy.l'),('U candidate_relative_to_reference/1','relativeEnergy.h')])
    s=s.replace('std::cout<<std::setprecision(24)',dump+'std::cout<<std::setprecision(24)')
   elif label.endswith('barrier'):
    dump=sidecar([('L marked_lower','lower'),('U longitudinal_upper','upper.h'),('L lateral_lower','lateral'),('RAW compiled_rho','rho'),('RAW compiled_delta','delta'),('RAW compiled_mustar','mustar')])
    s=s.replace('std::cout<<std::setprecision(24)',dump+'std::cout<<std::setprecision(24)')
   elif file=='certify_join_points.cpp':
    s=s.replace('std::cout<<std::setprecision(24)', 'std::ofstream ee("results/join_exact.hex");ee<<std::hexfloat<<"I q_gram "<<G.l<<" "<<G.h<<"\\n";\n std::cout<<std::setprecision(24)')
    anchor='if(k)std::cout<<",";std::cout<<'
    dump='ee<<"I points/"<<k<<"/distance_to_saddle_candidate "<<distance.l<<" "<<distance.h<<"\\nI points/"<<k<<"/signed_projection_on_normalized_q "<<p.l<<" "<<p.h<<"\\nU points/"<<k<<"/transverse_norm_upper "<<transverseUpper<<"\\nI points/"<<k<<"/relative_energy_to_reference "<<energy[k].l<<" "<<energy[k].h<<"\\n";\n '
    need(anchor in s,'join insertion missing');s=s.replace(anchor,dump+anchor)
    s=s.replace('return 0;','ee.close();need((bool)ee,"join exact sidecar output");report_guards();return 0;')
   elif file=='certify_path.cpp':
    s=s.replace('std::ofstream nodeTable(', 'std::ofstream ee(table+".exact.hex");ee<<std::hexfloat;\n std::ofstream nodeTable(')
    anchor='maximumGradientNorm=std::max(maximumGradientNorm,normup(norm));'
    s=s.replace(anchor,'ee<<"NODE "<<k<<" "<<energy.l<<" "<<energy.h<<" "<<normup(norm)<<"\\n";\n '+anchor)
    anchor='segmentTable<<"{'
    s=s.replace(anchor,'ee<<"SEG "<<k-1<<" "<<curv<<" "<<upper<<"\\n";'+anchor)
    s=s.replace('std::cout<<std::setprecision(24)', 'ee<<"U path_relative_reference_upper "<<maximum<<"\\n";\n std::cout<<std::setprecision(24)')
    s=s.replace('return 0;','ee.close();need((bool)ee,"path exact sidecar output");report_guards();return 0;')
   if label.startswith('stage2/') and not label.endswith('baseline'):s=s.replace('return 0;','report_guards();return 0;')
   (dst/file).write_text(s);(dst/(file+'.patch')).write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='original/'+label+'/'+file,tofile='instrumented/'+label+'/'+file)))
  subprocess.run(['g++',*FLAGS,file,'-o',exe],cwd=dst,check=True)
  records.append({'program':label+'/'+exe,'source':file,'flags':FLAGS})
 # Hardlink immutable proposals, never overwrite/write them.
 for p in src.rglob('*'):
  if p.is_file() and p.suffix not in ('.cpp',) and p.name not in ('MANIFEST.sha256',):
   rel=p.relative_to(src);target=dst/rel;target.parent.mkdir(parents=True,exist_ok=True)
   if not target.exists():target.hardlink_to(p)
 (dst/'results').mkdir(exist_ok=True)
base=SRC/'stage2/common/base_h.raw'
for target in (ROOT/'stage2/baseline/min.h',ROOT/'stage2/saddle/inputs/reference_h.raw',ROOT/'stage2/barrier/inputs/h.raw'):
 if not target.exists():target.hardlink_to(base)
for name in ('leg_return','leg_annihilation','full_annihilation'):
 p=ROOT/'connections'/(name+'.reference_h.raw')
 if not p.exists():p.hardlink_to(base)
N=256;nodebytes=N*N*8
with (ROOT/'connections/leg_return.z.raw').open('rb') as f:
 first=f.read(nodebytes);need(first==bytes(nodebytes),'return leg initial point differs');f.seek(-nodebytes,2);plus=f.read(nodebytes)
with (ROOT/'connections/leg_annihilation.z.raw').open('rb') as f:minus=f.read(nodebytes)
(ROOT/'connections/c_plus.z.raw').write_bytes(plus);(ROOT/'connections/c_minus.z.raw').write_bytes(minus)
(ROOT/'instrumentation.json').write_text(json.dumps(records,indent=2)+'\n')
print('CONNECTION_INSTRUMENTED_COPIES_BUILT')
