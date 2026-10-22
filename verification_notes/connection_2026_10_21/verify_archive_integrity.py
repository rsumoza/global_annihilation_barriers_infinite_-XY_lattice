"""Recheck every original member against FILES.sha256 and pinned ZIP identity."""
import argparse,hashlib,json,zipfile
from pathlib import Path
PINS=[('XY_Separacion_Global_Reproducible_2026-10-19.zip',228006035,'47eaf4886bc10b7bbb7ef259c44bfef9134bf3402920594ca6e20cd87235cedc'),('XY_Conexion_Global_Reproducible_2026-10-19.zip',441730856,'3a4163d8b24c21f95acbdf502ceaed46e967180c21abb363d3436ae5e2afb042'),('XY_Optimalidad_Extension_Reproducible_2026-10-19.zip',59642630,'42034c997678fa511ff261667ccb18edceee6d6bb05e543ddb6604f82958615d')]
def check(ok,why):
 if not ok: raise ValueError(why)
def main():
 p=argparse.ArgumentParser();p.add_argument('--archives',nargs=3,required=True);p.add_argument('--output',required=True);args=p.parse_args(); rows=[]
 for source,(name,size,want) in zip(args.archives,PINS):
  source=Path(source); actual=hashlib.sha256(source.read_bytes()).hexdigest();check(source.stat().st_size==size and actual==want,'Archive identity mismatch: '+name)
  with zipfile.ZipFile(source) as z:
   check(z.testzip() is None,'ZIP CRC failure: '+name); files=[n for n in z.namelist() if not n.endswith('/')]; manifests=[n for n in files if n.endswith('/FILES.sha256')];check(len(manifests)==1,'Expected unique FILES.sha256'); manifest=manifests[0]; root=manifest.rsplit('/',1)[0]+'/'; seen=[]
   for line in z.read(manifest).decode().splitlines():
    if not line.strip():continue
    digest,rel=line.split(maxsplit=1);rel=rel.lstrip('*'); entry=root+rel;check(entry in files,'Missing manifested member '+entry);check(hashlib.sha256(z.read(entry)).hexdigest()==digest,'Member hash mismatch '+entry);seen.append(entry)
   check(len(seen)==len(set(seen)),'Duplicate manifest member');check(set(files)==set(seen)|{manifest},'Unmanifested original member'); rows.append({'archive':name,'bytes':size,'sha256':actual,'members':len(files),'manifest_members':len(seen),'mismatches':[]})
 result={'status':'ORIGINAL_ARCHIVES_ALL_MEMBERS_PASS','archives':rows,'members_checked':sum(r['members'] for r in rows),'manifest_members_checked':sum(r['manifest_members'] for r in rows),'original_members_modified':False};Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
