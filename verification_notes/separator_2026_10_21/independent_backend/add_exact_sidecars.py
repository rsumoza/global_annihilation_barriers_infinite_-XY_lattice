from pathlib import Path
import subprocess,os,difflib,json,hashlib
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'xy_global_separator'/'distribution'
flags=['-std=c++17','-O2','-DNDEBUG','-fno-fast-math','-ffp-contract=off','-frounding-math']
s=(ROOT/'chart/certify.cpp').read_text()
new=r'''std::ofstream ee("results/coefficients_exact.hex");ee<<std::hexfloat;
 ee<<"HEADER mu "<<mu<<"\nHEADER delta "<<delta<<"\nHEADER remE "<<remE<<"\nHEADER remF "<<remF<<"\n";
 for(int k=0;k<2;k++)ee<<"A "<<k<<" "<<ac[k].l<<" "<<ac[k].h<<"\n";
 ee<<"Y2 "<<ynorm2.l<<" "<<ynorm2.h<<"\n";
 for(char tag:{'E','B'}){const Poly&pp=tag=='E'?ep:bp;for(int i=0;i<9;i++)for(int j=0;j<9-i;j++)if(pp[i][j].l||pp[i][j].h)ee<<tag<<" "<<i<<" "<<j<<" "<<pp[i][j].l<<" "<<pp[i][j].h<<"\n";}ee.close();need((bool)ee,"exact coefficient sidecar output");
 '''
s=s.replace('std::ofstream co("results/coefficients.txt");',new+'std::ofstream co("results/coefficients.txt");')
(ROOT/'chart/certify.cpp').write_text(s)
s=(ROOT/'geometry/certify.cpp').read_text()
new=r'''std::ofstream ee("geometry_exact.hex");ee<<std::hexfloat;
 for(int k=0;k<4;k++)ee<<"G "<<k<<" "<<G[k].l<<" "<<G[k].h<<"\n";
 for(int k=0;k<4;k++)ee<<"Gi "<<k<<" "<<Gi[k].l<<" "<<Gi[k].h<<"\n";
 ee<<"MAX neighbor "<<maxDiff<<"\nMAX Y "<<maxY<<"\nMAX coord "<<maxCoordinateError<<"\n";
 for(int k=0;k<count;k++){auto&c=C[k];for(int a=0;a<2;a++)ee<<"A "<<k<<" "<<a<<" "<<c.a[a].l<<" "<<c.a[a].h<<"\n";ee<<"Y2 "<<k<<" "<<c.ynorm2.l<<" "<<c.ynorm2.h<<"\n";}
 for(auto[i,j,d]:distances)ee<<"D "<<i<<" "<<j<<" "<<d.l<<" "<<d.h<<"\n";
 ee.close();need((bool)ee,"exact geometry sidecar output");
 '''
s=s.replace('std::cout<<std::setprecision(30)',new+'std::cout<<std::setprecision(30)',1)
(ROOT/'geometry/certify.cpp').write_text(s)
for name in ('chart','geometry'):
 d=ROOT/name
 subprocess.run(['g++',*flags,'certify.cpp','-o','verify_guarded_exact'],cwd=d,check=True)
 os.replace(d/'verify_guarded_exact',d/'verify_guarded')
 original=(SRC/name/'certify.cpp').read_text()
 (d/'certify.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),(d/'certify.cpp').read_text().splitlines(True),fromfile=f'original/{name}/certify.cpp',tofile=f'instrumented/{name}/certify.cpp')))
print('exact-sidecar sources atomically installed; scientific stdout unchanged')
