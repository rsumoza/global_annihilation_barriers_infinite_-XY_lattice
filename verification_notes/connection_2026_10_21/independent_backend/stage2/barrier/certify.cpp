#define main original_chain_main
#include "engine.cpp"
#undef main
using Poly=std::array<std::array<I,9>,9>;
static I eval(const Poly&p,I x,I y){std::array<I,9>xp,yp;xp[0]=yp[0]=I(1);for(int i=1;i<9;i++){xp[i]=xp[i-1]*x;yp[i]=yp[i-1]*y;}I s(0);for(int i=0;i<9;i++)for(int j=0;j<9-i;j++)if(p[i][j].l!=0||p[i][j].h!=0)s=s+p[i][j]*xp[i]*yp[j];return s;}
static int fac(int n){int a=1;for(int i=2;i<=n;i++)a*=i;return a;}
int main(int argc,char**argv){try{
 need(LDBL_MANT_DIG>=64&&std::numeric_limits<double>::is_iec559&&sizeof(double)==8,"platform");need(std::fegetround()==FE_TONEAREST,"rounding");uint16_t en=1;need(*(unsigned char*)&en==1,"endian");
 std::ifstream cfg("inputs/cfg");std::string st,rt,at,kt;cfg>>st>>rt>>at>>kt;need((bool)cfg,"cfg");R S=std::stold(st),Rr=std::stold(rt),a=std::stold(at),kap=std::stold(kt);need(S>.03L&&S<.6L&&Rr>.03L&&Rr<.4L&&a>0&&a<1&&kap>0&&kap<1,"domain");const R rho=.017422L,delta=.001166548L;
 auto h=readraw("inputs/h.raw",2048ull*2048),q=readraw("inputs/q.raw",8192);R hf=0;
 for(int i=0;i<2048;i++)for(int j=0;j<2048;j++)if(!(i>=993&&i<1055&&j>=993&&j<1055))hf=std::max(hf,std::fabs((R)h[(size_t)i*2048+j]));
 for(int k=0;k<2;k++)for(int i=0;i<64;i++)for(int j=0;j<64;j++)if(i==0||j==0||i==63||j==63)need(q[k*4096+i*64+j]==0,"q support");
 std::vector<CS>base(4096);for(int i=0;i<64;i++)for(int j=0;j<64;j++)base[i*64+j]=spin(i-31.5L+6,j-31.5L,12,h[(size_t)(i+992)*2048+j+992]);
 std::array<std::vector<I>,2>dq;for(auto&v:dq)v.resize(2*4096);std::vector<CS>bond;std::vector<Edge>all;
 I G00(0),G01(0),G11(0),H00(0),H01(0),H11(0);R maxLong=0,maxWidth=0;
 for(int ax=0;ax<2;ax++)for(int i=0;i<64-(ax==0);i++)for(int j=0;j<64-(ax==1);j++){
  int v=i*64+j,t=(i+(ax==0))*64+j+(ax==1),idx=ax*4096+v;
  dq[0][idx]=I(q[t])-I(q[v]);dq[1][idx]=I(q[4096+t])-I(q[4096+v]);G00=G00+sq(dq[0][idx]);G01=G01+dq[0][idx]*dq[1][idx];G11=G11+sq(dq[1][idx]);
  R lw=(I(S)*(I(mag(dq[0][idx]))+I(mag(dq[1][idx])))).h;maxLong=std::max(maxLong,lw);maxWidth=std::max(maxWidth,(I(lw)+I(Rr)/root(I(2))).h);
  CS bz=product(base[t],base[v]);H00=H00+bz.c*sq(dq[0][idx]);H01=H01+bz.c*dq[0][idx]*dq[1][idx];H11=H11+bz.c*sq(dq[1][idx]);bond.push_back(bz);all.push_back({ax,i,j});
 }
 R grammin=std::min((G00-I(mag(G01))).l,(G11-I(mag(G01))).l);need(grammin>0,"Gram lower");R rayleigh=(I(std::max((H00+I(mag(H01))).h,(H11+I(mag(H01))).h))/I(grammin)).h;need(rayleigh<.085L,"reaction subspace not soft");
 I det=G00*G11-sq(G01);need(G00.l>0&&det.l>0,"coordinate independence");std::array<I,4>Gi{G11/det,-G01/det,-G01/det,G00/det};
 need((root(Gi[0])*I(rho)).h<S&&(root(Gi[3])*I(rho)).h<S&&rho<Rr,"initial equilibrium outside reaction domain");need(maxWidth<1.12L,"domain does not preserve charges");
 I outang=I(12)/I(576)+I(2)*I(hf)+I(Rr)/root(I(2));need(outang.h<.4L&&trig(outang).c.l>a,"outer transverse Hessian");
 std::vector<Edge>ed;std::vector<I>weights;std::vector<std::array<I,2>>badB;
 for(size_t k=0;k<all.size();k++){
  Edge e=all[k];int idx=e.axis*4096+e.i*64+e.j;R lw=(I(S)*(I(mag(dq[0][idx]))+I(mag(dq[1][idx])))).h;I wid=I(lw)+I(Rr)/root(I(2));CS tr=trig(wid);CS z=bond[k];I lower=I(z.c.l)*(z.c.l>=0?tr.c:I(1))-I(mag(z.s))*tr.s;R lo=std::max((R)-1,lower.l);
  if(lo<a){ed.push_back(e);weights.push_back(I(a)-I(lo));badB.push_back({dq[0][idx],dq[1][idx]});}
 }
 std::ifstream ef("inputs/edges");need((bool)ef,"edge input open");int m=0;need((bool)(ef>>m),"edge count read");need(m==(int)ed.size()&&m>0&&m<300,"edge count");for(Edge e:ed){int ax=0,i=0,j=0;need((bool)(ef>>ax>>i>>j),"edge row read");need(ax==e.axis&&i==e.i&&j==e.j,"edge mismatch");}need((bool)ef,"edges");
 auto u=readraw("inputs/green.raw",(size_t)m*4096),V=readraw("inputs/basis.raw",(size_t)m*m);std::vector<R>err(m);std::vector<I>Rm(m*m),M(m*m);
 for(int k=0;k<m;k++){R s=0;auto[hi,hj]=head(ed[k]);for(int i=0;i<64;i++)for(int j=0;j<64;j++){int charge=(i==hi&&j==hj)-(i==ed[k].i&&j==ed[k].j);s=normsq_up_add(s,I(charge)-lap_at(&u[k*4096],64,i,j));}err[k]=(I(1024)*I(normup(s))).h;need(err[k]<1e-6L,"Green error");}
 for(int i=0;i<m;i++)for(int j=0;j<m;j++){auto[x,y]=head(ed[i]);R er=(root(I(2))*I(err[j])).h;Rm[i*m+j]=I(u[(size_t)j*4096+x*64+y])-I(u[(size_t)j*4096+ed[i].i*64+ed[i].j])+I(-er,er);}
 for(int i=0;i<m;i++)for(int j=0;j<i;j++){I v=meet(Rm[i*m+j],Rm[j*m+i]);Rm[i*m+j]=Rm[j*m+i]=v;}
 for(int i=0;i<m;i++)for(int j=0;j<m;j++){I projection(0);for(int k=0;k<2;k++)for(int l=0;l<2;l++)projection=projection+badB[i][k]*Gi[k*2+l]*badB[j][l];M[i*m+j]=I(i==j?kap:0)-root(weights[i])*(Rm[i*m+j]-projection)*root(weights[j])/I(a);}
 R dominance=positive_by_congruence(M,V,m),mu=(I(a)*(I(1)-I(kap))).l;need(mu>0,"transverse coercivity");std::cerr<<"transverse coercivity certified "<<mu<<"; polynomial residual next\n";
 // Taylor coefficients of the exact finite-support energy (degree 6) and force difference (degree 4).
 std::vector<std::array<int,2>>multi;for(int n=1;n<=4;n++)for(int i=0;i<=n;i++)multi.push_back({i,n-i});int tcount=multi.size();std::vector<std::vector<I>>fc(tcount,std::vector<I>(4096));Poly ep{};
 for(size_t k=0;k<all.size();k++){
  Edge e=all[k];auto[hi,hj]=head(e);int v=e.i*64+e.j,t=hi*64+hj,idx=e.axis*4096+v;I d0=dq[0][idx],d1=dq[1][idx];CS z=bond[k];
  for(int n=1;n<=6;n++)for(int i=0;i<=n;i++){int j=n-i;I factor=pwr(d0,i)*pwr(d1,j)/I(fac(i)*fac(j));I c=n%2?z.s:I(-1)*z.c;int sign=((n%2?(n-1)/2:n/2)%2)?-1:1;ep[i][j]=ep[i][j]+I(sign)*c*factor;}
  for(int tnum=0;tnum<tcount;tnum++){int i=multi[tnum][0],j=multi[tnum][1],n=i+j;I c=n%2?z.c:z.s;int sign=((n%2?(n-1)/2:n/2)%2)?-1:1;I flux=I(sign)*c*pwr(d0,i)*pwr(d1,j)/I(fac(i)*fac(j));fc[tnum][v]=fc[tnum][v]-flux;fc[tnum][t]=fc[tnum][t]+flux;}
 }
 // Each coefficient charge has exact mean zero by its divergence construction.
 auto pu=readraw("inputs/poly_green.raw",(size_t)tcount*4096);std::vector<R>pe(tcount),fn(tcount);std::vector<std::array<I,2>>proj(tcount);std::vector<I>Gram(tcount*tcount);
 for(int k=0;k<tcount;k++){R res=0,ns=0;for(int i=0;i<64;i++)for(int j=0;j<64;j++){int v=i*64+j;res=normsq_up_add(res,fc[k][v]-lap_at(&pu[k*4096],64,i,j));ns=normsq_up_add(ns,fc[k][v]);for(int l=0;l<2;l++)proj[k][l]=proj[k][l]+I(q[l*4096+v])*fc[k][v];}pe[k]=(I(1024)*I(normup(res))).h;fn[k]=normup(ns);need(pe[k]<1e-6L,"polynomial Poisson proposal inaccurate");}
 for(int i=0;i<tcount;i++)for(int j=0;j<tcount;j++){I v(0);for(int k=0;k<4096;k++)v=v+fc[i][k]*I(pu[j*4096+k]);R er=(I(fn[i])*I(pe[j])).h;Gram[i*tcount+j]=v+I(-er,er);}
 for(int i=0;i<tcount;i++)for(int j=0;j<i;j++){I v=meet(Gram[i*tcount+j],Gram[j*tcount+i]);Gram[i*tcount+j]=Gram[j*tcount+i]=v;}
 Poly bp{};for(int i=0;i<tcount;i++)for(int j=0;j<tcount;j++){I g=Gram[i*tcount+j];for(int k=0;k<2;k++)for(int l=0;l<2;l++)g=g-proj[i][k]*Gi[k*2+l]*proj[j][l];bp[multi[i][0]+multi[j][0]][multi[i][1]+multi[j][1]]=bp[multi[i][0]+multi[j][0]][multi[i][1]+multi[j][1]]+g;}
 // Uniform remainder bounds in the whole coordinate box.
 R gramnorm=std::max((G00+I(mag(G01))).h,(G11+I(mag(G01))).h);I qnorm=root(I(2)*sq(I(S))*I(gramnorm));R remF=(pwr(I(maxLong),4)*qnorm/I(120)).h,remE=(pwr(I(maxLong),5)*sq(qnorm)/I(5040)).h;
 auto bounds=[&](I x,I y){I ve=eval(ep,x,y)-I(remE);I b2=eval(bp,x,y);need(b2.h>=0,"negative upper bound for dual square");R beta=(I(delta)+root(I(std::max(std::ldexp(1.0L,-120),b2.h)))+I(remF)).h;return std::pair<I,R>{ve,beta};};
 // A-posteriori energy offset: Ehat-E* <= delta^2/(2*.05825).
 const R mustar=.05825L;I gap=sq(I(delta))/(I(2)*I(mustar));need(gap.h<.0000116811L,"energy gap");
 auto deriv=[](const Poly&p,int x,int y){Poly d{};for(int i=x;i<9;i++)for(int j=y;j<9-i;j++){int f=1;for(int k=0;k<x;k++)f*=i-k;for(int k=0;k<y;k++)f*=j-k;d[i-x][j-y]=I(f)*p[i][j];}return d;};
 const Poly ex=deriv(ep,1,0),ey=deriv(ep,0,1),bx=deriv(bp,1,0),by=deriv(bp,0,1);
 auto tight=[&](const Poly&p,const Poly&dx,const Poly&dy,I x,I y){R mx=(x.l+x.h)/2,my=(y.l+y.h)/2;I centered=eval(p,I(mx),I(my))+(x-I(mx))*eval(dx,x,y)+(y-I(my))*eval(dy,x,y);return meet(eval(p,x,y),centered);};
 R betamax=0,lateral=1e10L;
 auto enclosure=[&](I x,I y){I ve=tight(ep,ex,ey,x,y);I b2=tight(bp,bx,by,x,y);need(b2.h>=0,"dual square upper invalid");I beta=I(delta)+root(I(std::max(std::ldexp(1.0L,-120),b2.h)))+I(remF);betamax=std::max(betamax,beta.h);lateral=std::min(lateral,(ve-I(remE)-I(beta.h)*I(Rr)+I(mu)*sq(I(Rr))/I(2)).l);return I((ve-I(remE)-sq(beta)/(I(2)*I(mu))).l,(ve+I(remE)).h);};
 // Coefficients are output only after all inherited coercivity/Poisson checks pass.
 std::ofstream coef("results/coefficients.txt");coef<<std::setprecision(28)<<S<<" "<<mu<<" "<<delta<<" "<<remE<<" "<<remF<<" "<<gap.h<<"\n";
 for(char tag:{'E','B'}){const Poly&p=tag=='E'?ep:bp;for(int i=0;i<9;i++)for(int j=0;j<9-i;j++)if(p[i][j].l!=0||p[i][j].h!=0)coef<<tag<<" "<<i<<" "<<j<<" "<<p[i][j].l<<" "<<p[i][j].h<<"\n";}
 coef.close();need((bool)coef,"coefficient output write");
 int n=argc>1?std::stoi(argv[1]):160;need(n>=16&&n<=2048,"grid size");
 const size_t count=(size_t)n*n;std::vector<R>lo(count),hi(count),face(count,1e10L),ld(count,1e10L),ud(count,1e10L);std::vector<int>prev(count,-1);R section=1e10L;
 R xstart=(root(Gi[0])*I(rho)).h,ystart=(root(Gi[3])*I(rho)).h;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){I x=I(-S)+I(2)*I(S)*I(i,i+1)/I(n),y=I(-S)+I(2)*I(S)*I(j,j+1)/I(n);I en=enclosure(x,y);lo[i*n+j]=en.l;hi[i*n+j]=en.h;if(i==0)face[i*n+j]=std::min(face[i*n+j],enclosure(I(-S),y).l);if(i+1==n)face[i*n+j]=std::min(face[i*n+j],enclosure(I(S),y).l);if(j==0)face[i*n+j]=std::min(face[i*n+j],enclosure(x,I(-S)).l);if(j+1==n)face[i*n+j]=std::min(face[i*n+j],enclosure(x,I(S)).l);section=std::min(section,face[i*n+j]);}
 need(betamax<(I(mu)*I(Rr)).l,"transverse minimizer not certified inside tube");
 // Lower: ANY actual a* lies in this certified initial rectangle. The union
 // graph uses closed-box intersections, including corners, hence overcovers all
 // continuum paths; its bottleneck is a safe lower bound for escape.
 using Node=std::pair<R,int>;std::priority_queue<Node,std::vector<Node>,std::greater<Node>>heap;
 int starts=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++){I x=I(-S)+I(2)*I(S)*I(i,i+1)/I(n),y=I(-S)+I(2)*I(S)*I(j,j+1)/I(n);if(x.l<=xstart&&x.h>=-xstart&&y.l<=ystart&&y.h>=-ystart){int k=i*n+j;ld[k]=lo[k];heap.push({ld[k],k});starts++;}}
 R lower=1e10L;int lk=-1;while(!heap.empty()){auto[d,k]=heap.top();heap.pop();if(d!=ld[k])continue;int i=k/n,j=k%n;if(d>=lower)break;if(i==0||j==0||i+1==n||j+1==n){R candidate=std::max(d,face[k]);if(candidate<lower){lower=candidate;lk=k;}}for(int di=-1;di<=1;di++)for(int dj=-1;dj<=1;dj++){int ii=i+di,jj=j+dj;if(ii<0||jj<0||ii>=n||jj>=n||(di==0&&dj==0))continue;int kk=ii*n+jj;R dd=std::max(d,lo[kk]);if(dd<ld[kk]){ld[kk]=dd;heap.push({dd,kk});}}}
 // Upper: genuine polygonal path through y=0 configurations, starting at a=0.
 // The bridge theta*->thetahat costs at most gap. Every accepted box contains
 // only points with V <= hi, so straight segments through its centers/faces stay
 // below max hi. Four-neighbor adjacency makes this explicit.
 while(!heap.empty())heap.pop();for(int i=0;i<n;i++)for(int j=0;j<n;j++){I x=I(-S)+I(2)*I(S)*I(i,i+1)/I(n),y=I(-S)+I(2)*I(S)*I(j,j+1)/I(n);if(x.l<=0&&x.h>=0&&y.l<=0&&y.h>=0){int k=i*n+j;ud[k]=std::max((R)0,hi[k]);heap.push({ud[k],k});}}
 R upperV=1e10L;int uk=-1;while(!heap.empty()){auto[d,k]=heap.top();heap.pop();if(d!=ud[k])continue;int i=k/n,j=k%n;if(i==0||j==0||i+1==n||j+1==n){upperV=d;uk=k;break;}for(auto[di,dj]:{std::pair<int,int>{1,0},{-1,0},{0,1},{0,-1}}){int ii=i+di,jj=j+dj;if(ii<0||jj<0||ii>=n||jj>=n)continue;int kk=ii*n+jj;R dd=std::max(d,hi[kk]);if(dd<ud[kk]){ud[kk]=dd;prev[kk]=k;heap.push({dd,kk});}}}
 // Both the outer and inner graph are required to reach boundary.
 need(lk>=0&&uk>=0&&lower<=upperV+gap.h,"graph consistency");I upper=I(upperV)+gap;need(lower>0,"positive barrier not resolved");need(upper.h<lateral,"lateral escape not excluded below route");
 std::vector<int>path;for(int k=uk;k!=-1;k=prev[k])path.push_back(k);std::reverse(path.begin(),path.end());std::ofstream pf("results/path_grid_"+std::to_string(n)+".txt");pf<<std::setprecision(28);for(int k:path)pf<<k/n<<" "<<k%n<<" "<<hi[k]<<"\n";
 std::ofstream ee("results/exact.hex");ee<<std::hexfloat;
 ee<<"L marked_lower "<<lower<<"\n";
 ee<<"U longitudinal_upper "<<upper.h<<"\n";
 ee<<"L lateral_lower "<<lateral<<"\n";
 ee<<"RAW compiled_rho "<<rho<<"\n";
 ee<<"RAW compiled_delta "<<delta<<"\n";
 ee<<"RAW compiled_mustar "<<mustar<<"\n";
 ee.close();need((bool)ee,"exact sidecar output");
 std::cout<<std::setprecision(24)<<"{\"status\":\"MARKED_SUBLEVEL_PASS\",\"mesh\":"<<n<<",\"cells\":"<<count<<",\"start_cells\":"<<starts<<",\"start_radius_x\":"<<xstart<<",\"start_radius_y\":"<<ystart<<",\"energy_offset_upper\":"<<gap.h<<",\"transverse_coercivity\":"<<mu<<",\"lateral_lower\":"<<lateral<<",\"transverse_force_upper\":"<<betamax<<",\"simple_boundary_lower\":"<<section<<",\"marked_lower\":"<<lower<<",\"longitudinal_upper\":"<<upper.h<<",\"route_cells\":"<<path.size()<<",\"lower_exit_cell\":["<<lk/n<<","<<lk%n<<"],\"upper_exit_cell\":["<<uk/n<<","<<uk%n<<"],\"scope\":\"first exit from certified reaction tube only\"}\n";
 report_guards();return 0;
 }catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
