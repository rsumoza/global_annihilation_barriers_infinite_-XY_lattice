#!/usr/bin/env python3
"""Small negative controls for the exact-export checker, using actual sidecars."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('extension_export_audit',ROOT/'check_exports.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
lines=(ROOT/'results/projection_plus.exact.hex').read_text().splitlines()
raw=next(line.split()[2:4] for line in lines if line.startswith('I 0 '))
a,b=map(module.hf,raw)
records=[]
for name,binary,printed in [
    ('lower_decimal_excludes_actual_binary_lower',raw,(a+F(1,10**12),b+1)),
    ('upper_decimal_excludes_actual_binary_upper',raw,(a-1,b-F(1,10**12))),
    ('reversed_binary_interval',list(reversed(raw)),(a-1,b+1)),
]:
    rejected=False;reason=''
    try:module.compare_interval(binary,printed,name)
    except RuntimeError as e:rejected=True;reason=str(e)
    if not rejected:raise RuntimeError('negative export accepted: '+name)
    records.append({'name':name,'rejected':True,'reason':reason})
report={'status':'EXACT_EXPORT_NEGATIVE_CONTROLS_PASS','scope':'Three deliberately invalid interval-interface cases; not a full arbitrary-input validation','actual_sidecar_raw_interval':raw,'controls':records}
(ROOT/'results/negative_controls.json').write_text(json.dumps(report,indent=2)+'\n')
normal_spec=importlib.util.spec_from_file_location('independent_normal',ROOT/'check_normal_exact.py')
normal=importlib.util.module_from_spec(normal_spec);normal_spec.loader.exec_module(normal)
normal_report=json.loads((ROOT/'results/exact_normal_report.json').read_text());row=normal_report['records'][0]
rejected=False;reason=''
try:normal.radical_enclosure(-F(row['g'][0]),F(row['projected_norm_squared']),['0','0'],'zero interval excluding actual normal')
except RuntimeError as e:rejected=True;reason=str(e)
if not rejected:raise RuntimeError('invalid normal enclosure accepted')
(ROOT/'results/exact_normal_negative_control.json').write_text(json.dumps({'status':'EXACT_NORMAL_NEGATIVE_CONTROL_PASS','deliberately_invalid_interval':['0','0'],'reason':reason},indent=2)+'\n')
print(report['status'])
