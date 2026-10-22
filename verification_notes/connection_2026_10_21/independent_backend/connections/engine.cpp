// Second implementation. No Python checker, FFT, optimizer, or eigenvalue is trusted.
// Outward intervals in IEEE extended precision; elementary sqrt proposals are proved
// by squared endpoint inequalities. Transcendentals use independently coded series.
#include <algorithm>
#include <array>
#include <cfloat>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using R=long double;
static void need(bool c,const std::string&m){if(!c)throw std::runtime_error(m);}
static unsigned long long guarded_scalar_calls=0,guarded_interval_constructions=0;
static R dn(R x){guarded_scalar_calls++;need(std::isfinite(x),"nonfinite dn input");R y=std::nextafter(x,-std::numeric_limits<R>::infinity());need(std::isfinite(y),"nonfinite dn output");return y;}
static R up(R x){guarded_scalar_calls++;need(std::isfinite(x),"nonfinite up input");R y=std::nextafter(x,std::numeric_limits<R>::infinity());need(std::isfinite(y),"nonfinite up output");return y;}
struct I{R l,h;I(R a=0):l(a),h(a){guarded_interval_constructions++;need(std::isfinite(a),"nonfinite scalar interval");} I(R a,R b):l(a),h(b){guarded_interval_constructions++;need(std::isfinite(a)&&std::isfinite(b)&&a<=b,"invalid interval endpoints");}};
static void report_guards(){std::cerr<<"GUARD_DIAGNOSTIC finite_interval_constructions="<<guarded_interval_constructions<<" finite_scalar_calls="<<guarded_scalar_calls<<"\n";}
static I operator+(I a,I b){return {dn(a.l+b.l),up(a.h+b.h)};}
static I operator-(I a){return {-a.h,-a.l};}
static I operator-(I a,I b){return a+(-b);}
static I operator*(I a,I b){
 if((a.l==0&&a.h==0)||(b.l==0&&b.h==0))return I(0);
 R v[4]={a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h};
 return {dn(*std::min_element(v,v+4)),up(*std::max_element(v,v+4))};
}
static I operator/(I a,I b){need(b.l>0||b.h<0,"division crosses zero");return a*I(dn(1/b.h),up(1/b.l));}
static R mag(I a){return std::max(std::fabs(a.l),std::fabs(a.h));}
static I meet(I a,I b){I c(std::max(a.l,b.l),std::min(a.h,b.h));need(c.l<=c.h,"empty intersection");return c;}
static I sq(I a){R hi=up(mag(a)*mag(a));if(a.l<=0&&a.h>=0)return {0,hi};R lo=std::min(std::fabs(a.l),std::fabs(a.h));return {dn(lo*lo),hi};}
static I root(I a){
 need(a.l>=0&&std::isfinite(a.h),"sqrt domain");
 auto lower=[](R x){if(x==0)return (R)0;R r=dn(std::sqrt(x));int k=0;while(up(r*r)>x){r=dn(r);need(++k<100,"sqrt lower proposal rejected");}need(r>=0,"negative root");return r;};
 auto upper=[](R x){if(x==0)return (R)0;R r=up(std::sqrt(x));int k=0;while(dn(r*r)<x){r=up(r);need(++k<100,"sqrt upper proposal rejected");}return r;};
 return {lower(a.l),upper(a.h)};
}
static I pwr(I a,int k){I r(1);for(int i=0;i<k;i++)r=r*a;return r;}
struct CS{I c,s;};
static CS small(I x){
 need(mag(x)<=0.25L,"small trig domain");if(x.l==0&&x.h==0)return {I(1),I(0)};
 I z=sq(x),ct(1),st=x,c=ct,s=st;
 for(int j=1;j<=8;j++){ct=-(ct*z)/I((2*j-1)*(2*j));st=-(st*z)/I((2*j)*(2*j+1));c=c+ct;s=s+st;}
 R e=std::ldexp(1.0L,-80);return {c+I(-e,e),s+I(-e,e)};
}
static CS trig(I x){
 need(mag(x)<8,"general trig domain");CS a=small(x/I(32));
 for(int k=0;k<5;k++){I s=I(2)*a.s*a.c;I c=I(1)-I(2)*sq(a.s);a={c,s};}return a;
}
static I atanI(I x){
 need(mag(x)<=1.001L,"atan ratio domain");I t=x/(I(1)+root(I(1)+sq(x)));need(mag(t)<0.42L,"atan reduction");I t2=sq(t),v=t,sum(0);
 for(int k=0;k<36;k++){sum=sum+v/I(2*k+1);v=-(v*t2);}R e=std::ldexp(1.0L,-90);return I(2)*(sum+I(-e,e));
}
static I PI(){return I(16)*atanI(I(1)/I(5))-I(4)*atanI(I(1)/I(239));}
static I angle(I x,I y,I pi){
 R ax=mag(x),ay=mag(y);
 if(ax>=ay){
  need(x.l>0||x.h<0,"ambiguous angle axis");I a=atanI(y/x);
  if(x.l>0)return a;
  if(y.l>0)return a+pi;
  need(y.h<0,"ambiguous negative axis");return a-pi;
 }
 need(y.l>0||y.h<0,"ambiguous angle ordinate");I a=atanI(x/y);
 return (y.l>0?pi/I(2):-pi/I(2))-a;
}
static I principal(CS a,I pi){
 if(a.c.l>0&&mag(a.s)<=0.1L){
  I z=sq(a.s);I poly=I(35)/I(1152);poly=I(5)/I(112)+z*poly;poly=I(3)/I(40)+z*poly;poly=I(1)/I(6)+z*poly;poly=I(1)+z*poly;
  R e=(I(2)*pwr(I(mag(a.s)),11)).h;return a.s*poly+I(-e,e);
 }
 return angle(a.c,a.s,pi);
}
static R normsq_up_add(R acc,I x){if(x.l==0&&x.h==0)return acc;return up(acc+up(mag(x)*mag(x)));}
static R normup(R sqsum){return root(I(sqsum)).h;}
static void scalar_ok(R x,const std::string&m){need(std::isfinite(x),m);}
static std::vector<double> readraw(const std::string&fn,size_t count){
 std::ifstream f(fn,std::ios::binary|std::ios::ate);need((bool)f,"missing "+fn);need(f.tellg()==(std::streamoff)(count*sizeof(double)),"length "+fn);f.seekg(0);std::vector<double>a(count);f.read((char*)a.data(),count*sizeof(double));need((bool)f,"read "+fn);for(double x:a)need(std::isfinite(x),"nonfinite "+fn);return a;
}
static I lap_at(const double*u,int n,int i,int j){I a(u[i*n+j]),s(0);if(i)s=s+(a-I(u[(i-1)*n+j]));if(j)s=s+(a-I(u[i*n+j-1]));if(i+1<n)s=s+(a-I(u[(i+1)*n+j]));if(j+1<n)s=s+(a-I(u[i*n+j+1]));return s;}
struct Edge{int axis,i,j;};
static std::pair<int,int> head(Edge e){return {e.i+(e.axis==0),e.j+(e.axis==1)};}
static CS product(CS b,CS a){return {b.c*a.c+b.s*a.s,b.s*a.c-b.c*a.s};}
static CS spin(R x,R y,int L,R h){
 I X(x),Y(y),Z(x-L),ra=root(sq(X)+sq(Y)),rb=root(sq(Z)+sq(Y)),den=ra*rb;
 CS p={(X*Z+sq(Y))/den,I(-L)*Y/den},rot=small(I(h));
 return {p.c*rot.c-p.s*rot.s,p.s*rot.c+p.c*rot.s};
}
static R positive_by_congruence(const std::vector<I>&A,const std::vector<double>&V,int m){
 need(A.size()==(size_t)m*m&&V.size()==(size_t)m*m,"congruence dimensions");
 for(const I& x:A)need(std::isfinite(x.l)&&std::isfinite(x.h)&&x.l<=x.h,"nonfinite congruence A");
 for(double x:V)need(std::isfinite(x),"nonfinite congruence V");
 std::vector<I> av(m*m),t(m*m);R invert=1e10L,pivot=1e10L;
 // Invertibility of V is proved by V^T V being strictly positive diagonally dominant.
 for(int i=0;i<m;i++){
  R rad=0;I diag(0);
  for(int j=0;j<m;j++){
   I s(0);for(int k=0;k<m;k++)s=s+I(V[k*m+i])*I(V[k*m+j]);
   if(i==j)diag=s;else rad=up(rad+mag(s));
  }
  invert=std::min(invert,dn(diag.l-rad));
 }
 need(invert>0,"basis invertibility not certified");
 for(int i=0;i<m;i++)for(int j=0;j<m;j++){I s(0);for(int k=0;k<m;k++)s=s+A[i*m+k]*I(V[k*m+j]);av[i*m+j]=s;}
 for(int i=0;i<m;i++)for(int j=0;j<m;j++){I s(0);for(int k=0;k<m;k++)s=s+I(V[k*m+i])*av[k*m+j];t[i*m+j]=s;}
 for(int i=0;i<m;i++){R r=0;for(int j=0;j<m;j++)if(i!=j)r=up(r+mag(t[i*m+j]));pivot=std::min(pivot,dn(t[i*m+i].l-r));}
 need(pivot>0,"congruence diagonal dominance not certified");return pivot;
}
struct ChargeResult{std::vector<std::array<int,3>> q;R minc;};
static ChargeResult charge_patch(const std::vector<CS>&p,int side,int shift,I pi){
 std::vector<I> b0((side-1)*side),b1(side*(side-1));R minc=1;
 for(int i=0;i<side;i++)for(int j=0;j<side;j++){
  if(i+1<side){CS w=product(p[(i+1)*side+j],p[i*side+j]);minc=std::min(minc,w.c.l);b0[i*side+j]=principal(w,pi);}
  if(j+1<side){CS w=product(p[i*side+j+1],p[i*side+j]);minc=std::min(minc,w.c.l);b1[i*(side-1)+j]=principal(w,pi);}
 }
 ChargeResult result{{},minc};
 for(int i=0;i+1<side;i++)for(int j=0;j+1<side;j++){
  I a=(b0[i*side+j]+b1[(i+1)*(side-1)+j]-b0[i*side+j+1]-b1[i*(side-1)+j])/(I(2)*pi);
  int q=(int)std::llround((a.l+a.h)/2);need(a.l>q-0.25L&&a.h<q+0.25L,"charge integer not isolated");if(q)result.q.push_back({i-(side-2)/2+shift,j-(side-2)/2,q});
 }return result;
}
struct ChainResult{R c,delta,rho,mu,lower,G,V,energylo,energyhi,oldpath,m0;int negative;};
static ChainResult run(const std::string&pre,const std::string&pathfile){
 std::ostringstream out;
 std::ifstream cfg(pre+".cfg");int L,n;std::string astr,kstr;cfg>>L>>n>>astr>>kstr;need((bool)cfg,"config");R a0=std::stold(astr),kap=std::stold(kstr);need(L>=4&&L<=32&&L%2==0&&n>=64&&n%2==0&&n<=4096&&n>2*((L+1)/2)+12&&a0>0&&a0<1&&kap>0&&kap<1,"config domain");
 const int m=n+4,b=n+2,N=64,off=(m-N)/2;I pi=PI();R margin=(L==12?1.12L:1.30L);auto h=readraw(pre+".h",(size_t)n*n);R H=0;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){R x=h[i*n+j];need(x==-(R)h[i*n+n-1-j],"corrector not exactly odd");need((i&&j&&i+1<n&&j+1<n)||x==0,"corrector boundary not zero");H=std::max(H,std::fabs(x));}
 need(H<0.2L,"corrector size outside series domain");need((pi/I(2)+I(1)/I(L/2.0L-0.5L)+I(2*H)).h<pi.l,"raw angular lift homotopy not proved");
 // Check all Green proposal residuals before the expensive global assembly.
 std::ifstream ef(pre+".edges");need((bool)ef,"edge input open");int count=0;need((bool)(ef>>count),"edge count read");need(count>0&&count<256,"edge count");std::vector<Edge> proposed(count);for(auto&e:proposed){need((bool)(ef>>e.axis>>e.i>>e.j),"edge row read");auto[a,j]=head(e);need(e.axis>=0&&e.axis<2&&e.i>=0&&e.j>=0&&a<64&&j<64,"edge endpoints");}need((bool)ef,"edge input");
 auto gu=readraw(pre+".green",(size_t)count*4096);auto basis=readraw(pre+".basis",(size_t)count*count);std::vector<R>errors(count);
 for(int k=0;k<count;k++){R sum=0;auto[hi,hj]=head(proposed[k]);for(int i=0;i<64;i++)for(int j=0;j<64;j++){int f=(i==hi&&j==hj)-(i==proposed[k].i&&j==proposed[k].j);I r=I(f)-lap_at(&gu[k*4096],64,i,j);sum=normsq_up_add(sum,r);}errors[k]=(I(1024)*I(normup(sum))).h;need(errors[k]<1e-6L,"Green proposal too inaccurate");}
 std::cerr<<"input/Green proposals checked; assembling "<<L<<" "<<n<<"\n";
 std::vector<CS> p((size_t)m*m);std::vector<I>F((size_t)b*b,I(0));
 auto coord=[&](int i){return (R)i-(m-1)/2.0L;};
 for(int i=0;i<m;i++)for(int j=0;j<m;j++){R hh=(i>=2&&j>=2&&i<n+2&&j<n+2)?h[(i-2)*n+j-2]:0;p[(size_t)i*m+j]=spin(coord(i)+L/2,coord(j),L,hh);}
 std::vector<Edge>edges;std::vector<I>ws;std::array<std::vector<I>,2> raw;raw[0].resize((size_t)(m-1)*m);raw[1].resize((size_t)m*(m-1));R minc=1;std::vector<R>negupper;
 auto force_add=[&](int i,int j,I v){if(i>=1&&j>=1&&i<m-1&&j<m-1){size_t k=(size_t)(i-1)*b+j-1;F[k]=F[k]+v;}};
 for(int ax=0;ax<2;ax++)for(int i=0;i<m-(ax==0);i++)for(int j=0;j<m-(ax==1);j++){
  int hi=i+(ax==0),hj=j+(ax==1);CS w=product(p[(size_t)hi*m+hj],p[(size_t)i*m+j]);minc=std::min(minc,w.c.l);if(w.c.h<0)negupper.push_back(w.c.h);
  if(w.c.l<a0){need(i>=off&&j>=off&&hi<off+64&&hj<off+64,"deficit outside fixed Gram box");edges.push_back({ax,i-off,j-off});ws.push_back(I(a0)-w.c);need(ws.back().l>0,"deficit threshold ambiguous");}
  force_add(i,j,-w.s);force_add(hi,hj,w.s);
  I beta=principal(w,pi);if(ax==1&&coord(j)==-0.5L&&coord(i)+L/2>0&&coord(i)+L/2<L)beta=beta-I(2)*pi;
  raw[ax][(size_t)i*(m-(ax==1))+j]=beta;
 }
 need((int)edges.size()==count,"proposed edge set incomplete");for(int k=0;k<count;k++)need(edges[k].axis==proposed[k].axis&&edges[k].i==proposed[k].i&&edges[k].j==proposed[k].j,"proposed edges mismatch");
 // Outer weights: maximum phase increment, no missing bonds.
 int a=(L+1)/2,S=n/2+2;I outang=I(L)/sq(I(S-a-2));need(trig(outang).c.l>a0,"outer positivity");need(minc>(-trig(I(margin)).c).h,"global branch margin");
 std::vector<I>Rmat(count*count),K(count*count),A(count*count);
 for(int i=0;i<count;i++)for(int j=0;j<count;j++){auto[x,y]=head(edges[i]);I v=I(gu[(size_t)j*4096+x*64+y])-I(gu[(size_t)j*4096+edges[i].i*64+edges[i].j]);R er=(root(I(2))*I(errors[j])).h;Rmat[i*count+j]=v+I(-er,er);}
 for(int i=0;i<count;i++)for(int j=0;j<i;j++){I v=meet(Rmat[i*count+j],Rmat[j*count+i]);Rmat[i*count+j]=Rmat[j*count+i]=v;}
 for(int i=0;i<count;i++)for(int j=0;j<count;j++){K[i*count+j]=root(ws[i])*Rmat[i*count+j]*root(ws[j])/I(a0);A[i*count+j]=I(i==j?kap:0)-K[i*count+j];}
 R dominance=positive_by_congruence(A,basis,count);R c=(I(a0)*(I(1)-I(kap))).l;
 std::cerr<<"signed coercivity independently checked; residual next\n";
 auto u=readraw(pre+".u",(size_t)b*b);R energy=0,res2=0;
 for(int i=0;i<b;i++)for(int j=0;j<b;j++){
  if(i+1<b)energy=normsq_up_add(energy,I(u[(i+1)*b+j])-I(u[i*b+j]));if(j+1<b)energy=normsq_up_add(energy,I(u[i*b+j+1])-I(u[i*b+j]));
  res2=normsq_up_add(res2,F[(size_t)i*b+j]-lap_at(u.data(),b,i,j));
 }
 R finite=(I(normup(energy))+I(b)/I(2)*I(normup(res2))).h;S=b/2+1;I q=I(S)/I(S-a-2),sm(S-1);
 I tail=I(64)/I(5)*root(I(3))*I(L)*pwr(q,5)/(sq(sm)*root(sm))+I(32)/I(9)*root(I(3))*I(L*L*L)*pwr(q,7)/(pwr(sm,4)*root(sm));
 R delta=(I(finite)+I(tail.h)).h;I disc=sq(I(c))-I(2)*I(delta);need(disc.l>0,"Newton discriminant");I rr=I(2)*I(delta)/(I(c)+root(disc));R rho=up(rr.h*(1+std::ldexp(1.0L,-30)));need(((I(delta)+sq(I(rho))/I(2))/I(c)).h<=rho&&rho<c,"Newton invariant ball");R mu=(I(c)-I(rho)).l;need(rho<margin,"charges do not survive correction");int negatives=0;for(R w:negupper)if((I(w)+I(rho)).h<0)negatives++;
 std::vector<CS>patch(4096);for(int i=0;i<64;i++)for(int j=0;j<64;j++)patch[i*64+j]=p[(i+off)*m+j+off];
 auto charges=charge_patch(patch,64,L/2,pi);need(charges.q.size()==2,"initial charge count");bool plus=false,minus=false;for(auto v:charges.q){plus|=v==std::array<int,3>{0,0,1};minus|=v==std::array<int,3>{L,0,-1};}need(plus&&minus,"initial charge locations");
 need((I(L)/sq(I(30-a))+I(2*H)).h<pi.l,"unverified topology outside patch");
 R rad=std::min((I(2)*I(mu)).l,((I(margin)-I(rho))/I(2)).l);R lower=(I(mu)*sq(I(rad))/I(2)-pwr(I(rad),3)/I(6)).l;
 // Independently reproduce the complete global contraction path and displacement bounds.
 std::array<I,4>mom{I(0),I(0),I(0),I(0)};std::vector<I>high;R G2=0,V2=0;
 for(int ax=0;ax<2;ax++)for(int i=2;i<n+2;i++)for(int j=2;j<n+2;j++){
  int width=m-(ax==1);I d=raw[ax][(size_t)i*width+j];G2=normsq_up_add(G2,d);V2=normsq_up_add(V2,d-raw[ax][(size_t)(i-1)*width+j]);
  if(mag(d)<0.05L){I s=sq(d),v=s;for(int k=0;k<4;k++){mom[k]=mom[k]+v;v=v*s;}}else high.push_back(d);
 }
 S=n/2+1;q=I(S)/I(S-a-2);I tailE=I(4*L*L)*pwr(q,4)/sq(I(S-1));G2=(I(G2)+I(2)*tailE).h;R G=normup(G2);I qv=I(S)/I(S-a-3);V2=(I(V2)+I(32*L*L)*pwr(qv,4)/sq(I(S-1))).h;R V=normup(V2);
 I initialE;R maxE=-1e10L;
 for(int j=0;j<=64;j++){
  I t=I(j)/I(64),t2=sq(t),t4=sq(t2),t6=t4*t2,t8=sq(t4);I E=t2*mom[0]/I(2)-t4*mom[1]/I(24)+t6*mom[2]/I(720);R er=(t8*I(std::max((R)0,mom[3].h))/I(40320)).h;E=E+I(-er,er);for(I d:high)E=E+I(1)-trig(t*d).c;if(j==64)initialE=E;maxE=std::max(maxE,E.h);
 }
 R oldpath=(I(maxE)-I(initialE.l)+tailE+I(G2)/I(8*64*64)+I(rho)*I(G)+I(delta)*I(rho)+sq(I(rho))).h;
 R disp=(sq(I(V)+I(2)*I(rho))/I(2)).h;
 out<<std::setprecision(22)<<"{\"L\":"<<L<<",\"n\":"<<n<<",\"long_double_mantissa\":"<<LDBL_MANT_DIG<<",\"bad_edges\":"<<count<<",\"congruence_margin\":"<<dominance<<",\"c\":"<<c<<",\"residual_finite\":"<<finite<<",\"residual_tail\":"<<tail.h<<",\"delta\":"<<delta<<",\"rho\":"<<rho<<",\"mu\":"<<mu<<",\"negative_edges_survive\":"<<negatives<<",\"charges\":[[0,0,1],["<<L<<",0,-1]],\"branch_margin\":"<<margin<<",\"energy_candidate\":["<<initialE.l<<","<<(initialE+tailE).h<<"],\"barrier_lower\":"<<lower<<",\"global_contraction_upper\":"<<oldpath<<",\"displacement_upper\":"<<disp;
 // Optional explicit local path, checked independently with all time intervals covered.
 if(!pathfile.empty()){
  need(L==12,"local path is for L12");std::ifstream pc(pathfile+".cfg");int nodes,side;pc>>nodes>>side;need((bool)pc&&nodes>=2&&nodes<10000&&side==64,"path dimensions");auto z=readraw(pathfile+".z",(size_t)nodes*4096);
  auto zv=[&](int k,int i,int j)->R {return i>=0&&j>=0&&i<64&&j<64?z[(size_t)k*4096+i*64+j]:0;};
  for(int k=0;k<nodes;k++)for(int i=0;i<64;i++)for(int j=0;j<64;j++)need((i&&j&&i<63&&j<63)||(zv(k,i,j)==0),"path boundary");for(int i=0;i<4096;i++)need(z[i]==0,"path initial point");
  std::array<std::vector<CS>,2>wb;wb[0].resize(63*64);wb[1].resize(64*63);for(int ax=0;ax<2;ax++)for(int i=0;i<64-(ax==0);i++)for(int j=0;j<64-(ax==1);j++)wb[ax][i*(64-(ax==1))+j]=product(patch[(i+(ax==0))*64+j+(ax==1)],patch[i*64+j]);
  R prevU=0,maxPath=0,maxCurv=0;int seg=0;
  for(int k=0;k<nodes;k++){
   I E(0);R R2=0,curv2=0;
   for(int ax=0;ax<2;ax++)for(int i=0;i<64-(ax==0);i++)for(int j=0;j<64-(ax==1);j++){
    int hi=i+(ax==0),hj=j+(ax==1);I d=I(zv(k,hi,hj))-I(zv(k,i,j));CS rot=trig(d),base=wb[ax][i*(64-(ax==1))+j];I one=I(1)-rot.c;
    E=E+base.c*one+base.s*rot.s;R2=up(R2+up(2*std::max((R)0,one.h)));
    if(k){I a=(I(zv(k,hi,hj))-I(zv(k-1,hi,hj)))-(I(zv(k,i,j))-I(zv(k-1,i,j)));curv2=normsq_up_add(curv2,a);}
   }
   R U=(E+I(rho)*I(normup(R2))).h;
   if(k){R curv=(I(curv2)/I(8)).h;maxCurv=std::max(maxCurv,curv);maxPath=std::max(maxPath,(I(std::max(U,prevU))+I(curv)).h);seg++;}
   prevU=U;
  }
  // Endpoint charges and an explicit compact integer lifting.
  std::vector<CS>end(4096);std::vector<I>theta(4096);for(int i=0;i<64;i++)for(int j=0;j<64;j++){
   R x=i-31.5L+L/2,y=j-31.5L;CS rot=trig(I(zv(nodes-1,i,j)));CS base=patch[i*64+j];end[i*64+j]={base.c*rot.c-base.s*rot.s,base.s*rot.c+base.c*rot.s};
   int gi=i+off,gj=j+off;R hh=h[(gi-2)*n+gj-2];theta[i*64+j]=angle(I(x),I(y),pi)-angle(I(x-L),I(y),pi)+I(hh)+I(zv(nodes-1,i,j));
  }
  auto cq=charge_patch(end,64,L/2,pi);need(cq.q.empty(),"endpoint still has charges");need(cq.minc>(-trig(I(margin)).c).h,"endpoint branch margin");
  std::array<std::vector<int>,2>slip;slip[0].resize(63*64);slip[1].resize(64*63);
  for(int ax=0;ax<2;ax++)for(int i=0;i<64-(ax==0);i++)for(int j=0;j<64-(ax==1);j++){
   int t=(i+(ax==0))*64+j+(ax==1),s=i*64+j;I rawd=theta[t]-theta[s];I be=principal(product(end[t],end[s]),pi);I turns=(rawd-be)/(I(2)*pi);int k=(int)std::llround((turns.l+turns.h)/2);need(turns.l>k-.25L&&turns.h<k+.25L,"endpoint slip not integer isolated");slip[ax][i*(64-(ax==1))+j]=k;
  }
  std::vector<int>integer(4096,0);for(int i=0;i<64;i++)for(int j=0;j<64;j++){
   if(j)integer[i*64+j]=integer[i*64+j-1]+slip[1][i*63+j-1];else if(i)integer[i*64]=integer[(i-1)*64]+slip[0][(i-1)*64];
  }
  int maxturn=0;for(int i=0;i<64;i++)for(int j=0;j<64;j++){
   int v=integer[i*64+j];maxturn=std::max(maxturn,std::abs(v));if(i==0||j==0||i==63||j==63)need(v==0,"integer lift not compact");
   if(i+1<64)need(integer[(i+1)*64+j]-v==slip[0][i*64+j],"integer lift horizontal mismatch");if(j+1<64)need(integer[i*64+j+1]-v==slip[1][i*63+j],"integer lift vertical mismatch");
  }
  need(maxPath<0.031205L,"published local path upper bound not reproduced");need(lower>0.0001317L,"published local barrier lower bound not reproduced");
  out<<",\"path_segments\":"<<seg<<",\"path_curvature_max\":"<<maxCurv<<",\"local_path_upper\":"<<maxPath<<",\"endpoint_charge_count\":0,\"compact_integer_lift_max\":"<<maxturn;
 }
 out<<",\"status\":\"INDEPENDENT_CHAIN_PASS\"}\n";std::cout<<out.str();return {c,delta,rho,mu,lower,G,V,initialE.l,(initialE+tailE).h,oldpath,margin,negatives};
}
static void selftest(){
 need(LDBL_MANT_DIG>=64,"extended precision required");need(std::fegetround()==FE_TONEAREST,"rounding mode");I t=I(-std::ldexp(1.0L,-54))-I(1);need(t.l<-1&&t.h<=-1,"inexact difference test");I r=root(I(2));need(up(r.l*r.l)<=2&&dn(r.h*r.h)>=2,"sqrt certificate");I p=PI();need(p.l>3.14159265358979L&&p.h<3.14159265358980L,"Machin interval");bool rejected=false;try{positive_by_congruence({I(-1)},{1.0},1);}catch(...){rejected=true;}need(rejected,"negative quadratic form accepted");std::cout<<"{\"selftest\":\"PASS\",\"mantissa\":"<<LDBL_MANT_DIG<<",\"pi_low\":"<<std::setprecision(24)<<p.l<<",\"pi_high\":"<<p.h<<"}\n";
}
static void oracle(){
 std::string op,al,ah,bl,bh;
 while(std::cin>>op>>al>>ah>>bl>>bh){
  I a(std::stold(al),std::stold(ah)),b(std::stold(bl),std::stold(bh));
  need(std::isfinite(a.l)&&std::isfinite(a.h)&&a.l<=a.h&&std::isfinite(b.l)&&std::isfinite(b.h)&&b.l<=b.h,"oracle input");
  I r;
  if(op=="add")r=a+b;else if(op=="sub")r=a-b;else if(op=="mul")r=a*b;else if(op=="div")r=a/b;else if(op=="sq")r=sq(a);else if(op=="sqrt")r=root(a);
  else if(op=="sin")r=trig(a).s;else if(op=="cos")r=trig(a).c;else if(op=="atan")r=atanI(a);else if(op=="pi")r=PI();else need(false,"oracle operation");
  std::cout<<std::hexfloat<<r.l<<" "<<r.h<<"\n";
 }
}
int main(int argc,char**argv){try{
 need(LDBL_MANT_DIG>=64&&std::numeric_limits<double>::is_iec559&&sizeof(double)==8,"platform arithmetic");need(std::fegetround()==FE_TONEAREST,"round-to-nearest required");uint16_t endian=1;need(*(unsigned char*)&endian==1,"little endian required");
 if(argc==2&&std::string(argv[1])=="--oracle"){oracle();return 0;}
 if(argc==2&&std::string(argv[1])=="--selftest"){selftest();return 0;}need(argc==2||argc==3,"usage: check_chain PREFIX [PATH_PREFIX]");run(argv[1],argc==3?argv[2]:"");return 0;
 }catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
