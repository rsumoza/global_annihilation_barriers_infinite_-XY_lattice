"""Same finite/order checks, lazy diagnostic messages; no numeric edit."""
from pathlib import Path
import re,difflib,subprocess,json,hashlib
R=Path(__file__).resolve().parent;source=R/'stage2/saddle';d=R/'stage2/saddle_fast';d.mkdir(exist_ok=True)
old=(source/'engine.cpp').read_text();s=old.replace('static unsigned long long guarded_scalar_calls=0,guarded_interval_constructions=0;', 'static void guard_need(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}\nstatic unsigned long long guarded_scalar_calls=0,guarded_interval_constructions=0;')
for message in ['nonfinite dn input','nonfinite dn output','nonfinite up input','nonfinite up output','nonfinite scalar interval','invalid interval endpoints','congruence dimensions','nonfinite congruence A','nonfinite congruence V']:
 pattern=r'need\(([^;]*?),"'+re.escape(message)+r'"\)';s,n=re.subn(pattern,lambda m:'guard_need('+m.group(1)+',"'+message+'")',s)
 if n!=1:raise RuntimeError('unexpected guard source schema '+message)
(d/'engine.cpp').write_text(s);(d/'engine_lazy_messages.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='guarded/engine.cpp',tofile='guarded_lazy_messages/engine.cpp')))
(d/'certify.cpp').write_bytes((source/'certify.cpp').read_bytes());(d/'inputs').mkdir(exist_ok=True);(d/'results').mkdir(exist_ok=True)
for p in (source/'inputs').iterdir():
 if p.is_file() and not (d/'inputs'/p.name).exists():(d/'inputs'/p.name).hardlink_to(p)
flags=json.loads((R/'instrumentation.json').read_text())[0]['flags']
subprocess.run(['g++',*flags,'certify.cpp','-o','verify_guarded_fast'],cwd=d,check=True)
(R/'lazy_message_source_binding.json').write_text(json.dumps({'scope':'only eager construction of diagnostic messages changed; predicates,numeric formulas,counters unchanged','original_guarded_engine_sha256':hashlib.sha256(old.encode()).hexdigest(),'lazy_message_guarded_engine_sha256':hashlib.sha256(s.encode()).hexdigest(),'certify_source_identical':(source/'certify.cpp').read_bytes()==(d/'certify.cpp').read_bytes(),'flags':flags},indent=2)+'\n')
print('LAZY_MESSAGE_GUARDED_SADDLE_BUILT')
