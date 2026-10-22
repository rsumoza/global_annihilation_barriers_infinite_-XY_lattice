#!/usr/bin/env python3
"""Reconstruct exactly the frozen extension inputs and compile guarded copies.

Original packages are only read. This is fixed-input diagnostic instrumentation,
not an independently implemented interval engine or an official source revision.
"""
from argparse import ArgumentParser
from pathlib import Path
import difflib
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parents[1]
FLAGS = ['-std=c++17', '-O2', '-DNDEBUG', '-fno-fast-math', '-ffp-contract=off', '-frounding-math']
PINS = {
    'extension': '248f0614524832610b07fd21b5d6b9eaf670d1846d6032aa22a51755a5c95972',
    'separator': '991c744eab7e812d629662a0c0bc5cbc177e9dd2e9daa67caeefab7a814f941b',
    'connection': 'ba100e9b7b8763bcd9ce1bcb7ee1a50a1f76669d2be92966ecbe6f6208feb27c',
}
CHANGED = (42, 43, 44, 53, 55)

def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)

def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def replace_once(s, old, new):
    need(s.count(old) == 1, 'source patch anchor count differs: ' + old[:80])
    return s.replace(old, new)

def ledger(root, pin):
    manifest = root / 'MANIFEST.sha256'
    need(sha(manifest) == pin, 'frozen manifest differs: ' + str(root))
    result = {}
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        p = Path(name)
        need(not p.is_absolute() and '..' not in p.parts and name not in result, 'unsafe/duplicate member')
        need(sha(root / p) == digest, 'source/data identity differs: ' + name)
        result[name] = digest
    return result

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--extension', type=Path, default=BASE/'recovered/unpacked/xy_optimality_extension')
    parser.add_argument('--separator', type=Path, default=BASE/'audit_separator_2026_10_21/xy_global_separator')
    parser.add_argument('--connection', type=Path, default=ROOT.parent/'xy_global_connection')
    args = parser.parse_args()
    extension = args.extension.resolve()/'extension'
    separator = args.separator.resolve()/'distribution'
    connection = args.connection.resolve()/'distribution'
    em, sm, cm = ledger(extension, PINS['extension']), ledger(separator, PINS['separator']), ledger(connection, PINS['connection'])
    work = ROOT/'work';work.mkdir(exist_ok=True)
    bindings = []
    def copy(src, root, name, dest, expected):
        need(sha(root/name) == expected, 'bound source changed: ' + name)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root/name, dest)
        need(sha(dest) == expected, 'bound copy differs')
        bindings.append({'source_package':src, 'path':name, 'sha256':expected, 'bytes':dest.stat().st_size})
    for name,digest in em.items():
        copy('extension',extension,name,work/name,digest)
    for root,book,name,dest in [
        (separator,sm,'common/h.raw','common/base_h.raw'),
        (separator,sm,'common/q.raw','common/Q_old.raw'),
        (connection,cm,'stage2/saddle/inputs/h.raw','common/saddle_h.raw'),
        (connection,cm,'stage2/saddle/inputs/q.raw','common/q_negative.raw'),
    ]:
        copy('separator' if root==separator else 'connection',root,name,work/dest,book[name])
    need(cm['stage2/common/base_h.raw']==sm['common/h.raw'], 'reference identity mismatch')
    for cid in range(109):
        name=f'centers/center_{cid:03}.raw';copy('separator',separator,name,work/name,sm[name])
    (ROOT/'input_bindings.json').write_text(json.dumps({'pins':PINS,'inputs':bindings},indent=2)+'\n')
    source = ROOT/'instrumented';source.mkdir(exist_ok=True)
    original = (extension/'engine.cpp').read_text();s = original
    s=replace_once(s,'static R dn(R x){return std::nextafter(x,-std::numeric_limits<R>::infinity());}',
'''static unsigned long long guarded_scalar_calls=0,guarded_interval_constructions=0;
static R dn(R x){guarded_scalar_calls++;need(std::isfinite(x),"nonfinite dn input");R y=std::nextafter(x,-std::numeric_limits<R>::infinity());need(std::isfinite(y),"nonfinite dn output");return y;}''')
    s=replace_once(s,'static R up(R x){return std::nextafter(x,std::numeric_limits<R>::infinity());}',
'''static R up(R x){guarded_scalar_calls++;need(std::isfinite(x),"nonfinite up input");R y=std::nextafter(x,std::numeric_limits<R>::infinity());need(std::isfinite(y),"nonfinite up output");return y;}''')
    s=replace_once(s,'struct I{R l,h;I(R a=0):l(a),h(a){} I(R a,R b):l(a),h(b){}};',
'''struct I{R l,h;I(R a=0):l(a),h(a){guarded_interval_constructions++;need(std::isfinite(a),"nonfinite scalar interval");} I(R a,R b):l(a),h(b){guarded_interval_constructions++;need(std::isfinite(a)&&std::isfinite(b)&&a<=b,"invalid interval endpoints");}};
static void report_guards(){std::cerr<<"GUARD_DIAGNOSTIC finite_interval_constructions="<<guarded_interval_constructions<<" finite_scalar_calls="<<guarded_scalar_calls<<"\\n";}
static std::ofstream exact_stream("exact_exports.hex");
static unsigned long long exact_record_count=0;
static void exact_interval(I a){exact_stream<<std::hexfloat<<"I "<<exact_record_count++<<" "<<a.l<<" "<<a.h<<"\\n";}
static void exact_constant(const char*name,R a){exact_stream<<std::hexfloat<<"CONST "<<name<<" "<<a<<"\\n";}
static void finish_exact(){exact_stream.flush();need((bool)exact_stream,"exact export sidecar failed");report_guards();}''')
    s=replace_once(s,'std::vector<I> av(m*m),t(m*m);R invert=1e10L,pivot=1e10L;',
'''need(A.size()==(size_t)m*m&&V.size()==(size_t)m*m,"congruence dimensions");
 for(const I& x:A)need(std::isfinite(x.l)&&std::isfinite(x.h)&&x.l<=x.h,"invalid congruence A");
 for(double x:V)need(std::isfinite(x),"nonfinite congruence V");
 std::vector<I> av(m*m),t(m*m);R invert=1e10L,pivot=1e10L;''')
    (source/'engine.cpp').write_text(s)
    (source/'engine.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile='original/engine.cpp',tofile='instrumented/engine.cpp')))
    sources = {'engine.cpp':{'original':sha(extension/'engine.cpp'),'guarded':sha(source/'engine.cpp')}}
    for name in ('projection','center_distances','adjoint','chart'):
        filename=f'certify_{name}.cpp';old=(extension/filename).read_text();s=old
        if name in ('projection','center_distances'):
            s=replace_once(s,'static void jout(I a){','static void jout(I a){exact_interval(a);')
        elif name=='adjoint':
            s=replace_once(s,'static void ji(I a){','static void ji(I a){exact_interval(a);')
        else:
            sidecar=r'''std::ofstream ee("results/coefficients_exact.hex");ee<<std::hexfloat;
 ee<<"HEADER mu "<<mu<<"\nHEADER delta "<<delta<<"\nHEADER remE "<<remE<<"\nHEADER remF "<<remF<<"\n";
 for(int k=0;k<2;k++)ee<<"A "<<k<<" "<<ac[k].l<<" "<<ac[k].h<<"\n";
 ee<<"Y2 "<<ynorm2.l<<" "<<ynorm2.h<<"\n";
 for(char tag:{'E','B'}){const Poly&pp=tag=='E'?ep:bp;for(int i=0;i<9;i++)for(int j=0;j<9-i;j++)if(pp[i][j].l||pp[i][j].h)ee<<tag<<" "<<i<<" "<<j<<" "<<pp[i][j].l<<" "<<pp[i][j].h<<"\n";}ee.close();need((bool)ee,"exact chart sidecar failed");
 '''
            s=replace_once(s,'std::ofstream co("results/coefficients.txt");',sidecar+'std::ofstream co("results/coefficients.txt");')
        if name=='projection':
            s=replace_once(s,'std::cout<<std::setprecision(30)',
                'exact_constant("rho_safe",.012825L);exact_constant("m_safe",.062886L);exact_constant("candidate_complement",.0719549999L);\n std::cout<<std::setprecision(30)')
        if name=='adjoint':
            s=replace_once(s,'std::cout<<std::setprecision(30)',
                'exact_constant("rho_safe",.012825L);\n std::cout<<std::setprecision(30)')
        s=replace_once(s,'return 0;','finish_exact();return 0;')
        (source/filename).write_text(s)
        (source/(filename+'.patch')).write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='original/'+filename,tofile='instrumented/'+filename)))
        with (ROOT/(name+'_build.stdout.log')).open('w') as out,(ROOT/(name+'_build.stderr.log')).open('w') as err:
            subprocess.run(['g++',*FLAGS,str(source/filename),'-o',str(work/('guarded_'+name))],stdout=out,stderr=err,check=True)
        sources[filename]={'original':sha(extension/filename),'guarded':sha(source/filename)}
    for cid in CHANGED:
        stem=f'c{cid:03}';inp=work/'chart_work'/stem/'inputs';inp.mkdir(parents=True,exist_ok=True);(inp.parent/'results').mkdir(exist_ok=True)
        for member in ('cfg','edges','basis.raw','green.raw','poly_green.raw'):
            shutil.copyfile(work/'cases'/stem/member,inp/member)
        for src,dest in [(work/'common/base_h.raw','h.raw'),(work/'common/Q_old.raw','q.raw'),(work/'centers'/f'center_{cid:03}.raw','center.raw')]:
            shutil.copyfile(src,inp/dest)
    (ROOT/'instrumentation.json').write_text(json.dumps({'flags':FLAGS,'source_hashes':sources,'scope':'Fixed extension inputs, guarded original arithmetic formulas plus exact export sidecars; originals unchanged'},indent=2)+'\n')
    print('Guarded extension reconstructed and compiled',flush=True)

if __name__=='__main__':main()
