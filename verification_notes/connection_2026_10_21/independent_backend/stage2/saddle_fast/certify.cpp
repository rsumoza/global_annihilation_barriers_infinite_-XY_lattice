#define main original_chain_main
#include "engine.cpp"
#undef main

// phi+h is represented by its unit complex phase. The full candidate has
// |h| around .408, so the old small-angle guard is not used for that rotation.
static CS candidate_spin(R x,R y,I h){
 I X(x),Y(y),Z(x-12),den=root(sq(X)+sq(Y))*root(sq(Z)+sq(Y));
 CS p{(X*Z+sq(Y))/den,I(-12)*Y/den};
 CS r=mag(h)<=.25L?small(h):trig(h);
 return {p.c*r.c-p.s*r.s,p.s*r.c+p.c*r.s};
}
int main(int argc,char**argv){try{
 need(argc==1||argc==2,"usage: certify [INPUT_DIRECTORY]");std::string input=argc==2?argv[1]:"inputs";
 need(LDBL_MANT_DIG>=64&&std::numeric_limits<double>::is_iec559&&sizeof(double)==8,"platform");
 need(std::fegetround()==FE_TONEAREST,"rounding");uint16_t endian=1;need(*(unsigned char*)&endian==1,"endian");
 std::ifstream cfg(input+"/cfg");int n=0;std::string as,ks;need((bool)(cfg>>n>>as>>ks),"config read");
 R a0=std::stold(as),kap=std::stold(ks);need((n==2048||n==4096)&&a0>0&&a0<1&&kap>0&&kap<1,"config domain");
 const int m=n+4,b=n+2,N=64,off=(m-N)/2;const R margin=.65L;I pi=PI();
 auto h=readraw(input+"/h.raw",(size_t)n*n),q=readraw(input+"/q.raw",4096);R H=0,Hfar=0;bool odd=true;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){
  R x=h[(size_t)i*n+j];H=std::max(H,std::fabs(x));odd&=x==-(R)h[(size_t)i*n+n-1-j];
  need((i&&j&&i+1<n&&j+1<n)||x==0,"candidate boundary not zero");
  if(!(i>=n/2-31&&i<n/2+31&&j>=n/2-31&&j<n/2+31))Hfar=std::max(Hfar,std::fabs(x));
 }
 need(H<.5L,"candidate phase domain");
 for(int i=0;i<64;i++)for(int j=0;j<64;j++)need((i&&j&&i<63&&j<63)||q[i*64+j]==0,"q support");
 std::array<std::vector<I>,2>dq;dq[0].resize(63*64);dq[1].resize(64*63);I G(0);
 for(int ax=0;ax<2;ax++)for(int i=0;i<64-(ax==0);i++)for(int j=0;j<64-(ax==1);j++){
  int width=64-(ax==1);I d=I(q[(i+(ax==0))*64+j+(ax==1)])-I(q[i*64+j]);dq[ax][i*width+j]=d;G=G+sq(d);
 }
 need(G.l>.5L&&G.h<2,"q nonzero and scale");
 std::ifstream ef(input+"/edges");need((bool)ef,"edge open");int count=0;need((bool)(ef>>count),"edge count read");need(count>0&&count<256,"edge count");std::vector<Edge>proposed(count);
 for(auto&e:proposed){need((bool)(ef>>e.axis>>e.i>>e.j),"edge row read");need(e.axis>=0&&e.axis<2&&e.i>=0&&e.i<64&&e.j>=0&&e.j<64,"edge range");auto[x,y]=head(e);need(x<64&&y<64,"edge endpoint");}
 auto gu=readraw(input+"/green.raw",(size_t)count*4096),V=readraw(input+"/basis.raw",(size_t)count*count);std::vector<R>err(count);
 for(int k=0;k<count;k++){
  R r2=0;auto[x,y]=head(proposed[k]);
  for(int i=0;i<64;i++)for(int j=0;j<64;j++)r2=normsq_up_add(r2,I((i==x&&j==y)-(i==proposed[k].i&&j==proposed[k].j))-lap_at(&gu[(size_t)k*4096],64,i,j));
  err[k]=(I(1024)*I(normup(r2))).h;need(err[k]<1e-6L,"Green proposal inaccurate");
 }
 std::cerr<<"Proposals checked; assembling the complete infinite-domain candidate\n";
 std::vector<CS>p((size_t)m*m);std::vector<I>F((size_t)b*b,I(0));
 auto coord=[&](int i){return (R)i-(m-1)/2.0L;};
 for(int i=0;i<m;i++)for(int j=0;j<m;j++){
  I hh((i>=2&&j>=2&&i<n+2&&j<n+2)?h[(size_t)(i-2)*n+j-2]:0);
  p[(size_t)i*m+j]=candidate_spin(coord(i)+6,coord(j),hh);
 }
 auto force_add=[&](int i,int j,I value){if(i>=1&&j>=1&&i<m-1&&j<m-1){size_t k=(size_t)(i-1)*b+j-1;F[k]=F[k]+value;}};
 auto reference=readraw(input+"/reference_h.raw",2048ull*2048);
 auto difference=[&](int i,int j){
  if(i<2||j<2||i>=n+2||j>=n+2)return I(0);int ri=i-2-(n-2048)/2,rj=j-2-(n-2048)/2;R old=(ri>=0&&rj>=0&&ri<2048&&rj<2048)?reference[(size_t)ri*2048+rj]:0;R now=h[(size_t)(i-2)*n+j-2];return now==old?I(0):I(now)-I(old);
 };
 std::vector<Edge>edges;std::vector<I>w,Bq;I Hq(0),E(0),relativeEnergy(0);R minc=1;
 for(int ax=0;ax<2;ax++)for(int i=0;i<m-(ax==0);i++)for(int j=0;j<m-(ax==1);j++){
  int hi=i+(ax==0),hj=j+(ax==1);CS z=product(p[(size_t)hi*m+hj],p[(size_t)i*m+j]);
  minc=std::min(minc,z.c.l);E=E+I(1)-z.c;force_add(i,j,-z.s);force_add(hi,hj,z.s);
  I dd=difference(hi,hj)-difference(i,j);CS dr=(dd.l==0&&dd.h==0)?CS{I(1),I(0)}:(mag(dd)<=.25L?small(dd):trig(dd));relativeEnergy=relativeEnergy+z.c*(dr.c-I(1))+z.s*dr.s;
  bool inside=i>=off&&j>=off&&hi<off+64&&hj<off+64;
  if(inside){I d=dq[ax][(i-off)*(64-(ax==1))+j-off];Hq=Hq+z.c*sq(d);}
  if(z.c.l<a0){need(inside,"deficient edge outside Gram box");edges.push_back({ax,i-off,j-off});w.push_back(I(a0)-z.c);need(w.back().l>0,"ambiguous deficient threshold");Bq.push_back(dq[ax][(i-off)*(64-(ax==1))+j-off]);}
 }
 need((int)edges.size()==count,"incomplete deficient edge list");for(int k=0;k<count;k++)need(edges[k].axis==proposed[k].axis&&edges[k].i==proposed[k].i&&edges[k].j==proposed[k].j,"deficient edge mismatch");
 int S=n/2+2;I outerAngle=I(12)/sq(I(S-6-2));need(trig(outerAngle).c.l>a0,"uncontrolled exterior Hessian");
 need(minc>(-trig(I(margin)).c).h,"candidate branch margin");
 std::vector<I>Rmat(count*count),A(count*count);
 for(int i=0;i<count;i++)for(int j=0;j<count;j++){
  auto[x,y]=head(edges[i]);I v=I(gu[(size_t)j*4096+x*64+y])-I(gu[(size_t)j*4096+edges[i].i*64+edges[i].j]);R error=(root(I(2))*I(err[j])).h;Rmat[i*count+j]=v+I(-error,error);
 }
 for(int i=0;i<count;i++)for(int j=0;j<i;j++){I v=meet(Rmat[i*count+j],Rmat[j*count+i]);Rmat[i*count+j]=Rmat[j*count+i]=v;}
 for(int i=0;i<count;i++)for(int j=0;j<count;j++)A[i*count+j]=I(i==j?kap:0)-root(w[i])*(Rmat[i*count+j]-Bq[i]*Bq[j]/G)*root(w[j])/I(a0);
 R dominance=positive_by_congruence(A,V,count),positive=(I(a0)*(I(1)-I(kap))).l;I alpha=Hq/G;need(alpha.h<0,"negative direction not certified");
 R negative=(-alpha).l,sigma=std::min(positive,negative);need(sigma>0,"indefinite invertibility");R inverse=(I(1)/I(sigma)).h;
 std::vector<CS>patch(4096);for(int i=0;i<64;i++)for(int j=0;j<64;j++)patch[i*64+j]=p[(size_t)(i+off)*m+j+off];
 auto charges=charge_patch(patch,64,6,pi);need(charges.q.size()==2,"charge count");bool plus=false,minus=false;for(auto v:charges.q){plus|=v==std::array<int,3>{0,0,1};minus|=v==std::array<int,3>{12,0,-1};}need(plus&&minus,"charge positions");
 need((I(12)/sq(I(24))+I(2)*I(Hfar)).h<pi.l,"uncontrolled exterior charges");
 std::cerr<<"Index one and exterior coercivity certified; validating global residual\n";
 auto u=readraw(input+"/u.raw",(size_t)b*b);R en=0,res=0;
 for(int i=0;i<b;i++)for(int j=0;j<b;j++){
  if(i+1<b)en=normsq_up_add(en,I(u[(size_t)(i+1)*b+j])-I(u[(size_t)i*b+j]));
  if(j+1<b)en=normsq_up_add(en,I(u[(size_t)i*b+j+1])-I(u[(size_t)i*b+j]));
  res=normsq_up_add(res,F[(size_t)i*b+j]-lap_at(u.data(),b,i,j));
 }
 // The exact finite force has mean zero: subtract the pure dipole force.
 // Its flux difference is compact because h vanishes on its boundary and
 // outside. The pure dipole force sums to zero by reflection in y.
 R finite=(I(normup(en))+I(b)/I(2)*I(normup(res))).h;S=b/2+1;I ratio=I(S)/I(S-8),sm(S-1);
 I tail=I(64)/I(5)*root(I(3))*I(12)*pwr(ratio,5)/(sq(sm)*root(sm))+I(32)/I(9)*root(I(3))*I(12*12*12)*pwr(ratio,7)/(pwr(sm,4)*root(sm));
 R delta=(I(finite)+tail).h;I disc=sq(I(sigma))-I(2)*I(delta);need(disc.l>0,"indefinite Newton discriminant");I sharp=I(2)*I(delta)/(I(sigma)+root(disc));R rho=up(sharp.h*(1+std::ldexp(1.0L,-30)));
 need(((I(delta)+sq(I(rho))/I(2))/I(sigma)).h<=rho&&rho<sigma,"indefinite Newton ball");need(rho<margin,"charges not preserved by correction");
 R exactPositive=(I(positive)-I(rho)).l,exactNegative=(I(negative)-I(rho)).l;need(exactPositive>0&&exactNegative>0,"index does not survive correction");
 S=n/2+1;I rat=I(S)/I(S-8),energyTail=I(4*12*12)*pwr(rat,4)/sq(I(S-1));
 I candidateEnergy(E.l,(E+energyTail).h);R energyError=std::min((sq(I(rho))/I(2)).h,(I(delta)*I(rho)/I(2)+pwr(I(rho),3)/I(12)).h);
 std::ofstream ee("results/exact.hex");ee<<std::hexfloat;
 ee<<"U delta "<<delta<<"\n";
 ee<<"L negative_direction_bound "<<negative<<"\n";
 ee<<"L transverse_coercivity "<<positive<<"\n";
 ee<<"U rho "<<rho<<"\n";
 ee<<"L candidate_branch_margin "<<margin<<"\n";
 ee<<"L candidate_relative_to_reference/0 "<<relativeEnergy.l<<"\n";
 ee<<"U candidate_relative_to_reference/1 "<<relativeEnergy.h<<"\n";
 ee.close();need((bool)ee,"exact sidecar output");
 std::cout<<std::setprecision(24)<<"{\"status\":\"INFINITE_SADDLE_INDEX_ONE_PASS\",\"n\":"<<n<<",\"q_support_side\":64,\"long_double_mantissa\":"<<LDBL_MANT_DIG<<",\"h_max\":"<<H<<",\"h_exact_y_odd_diagnostic\":"<<(odd?"true":"false")<<",\"mean_zero_proof\":\"compact_flux_plus_reflected_dipole\",\"deficient_edges\":"<<count<<",\"kappa\":"<<kap<<",\"congruence_margin\":"<<dominance<<",\"q_gram\":["<<G.l<<","<<G.h<<"],\"negative_rayleigh\":["<<alpha.l<<","<<alpha.h<<"],\"transverse_coercivity\":"<<positive<<",\"negative_direction_bound\":"<<negative<<",\"sigma\":"<<sigma<<",\"inverse_norm_upper\":"<<inverse<<",\"residual_finite\":"<<finite<<",\"residual_tail\":"<<tail.h<<",\"delta\":"<<delta<<",\"rho\":"<<rho<<",\"newton_contraction_upper\":"<<(I(rho)/I(sigma)).h<<",\"exact_transverse_coercivity\":"<<exactPositive<<",\"exact_negative_direction_bound\":"<<exactNegative<<",\"candidate_branch_margin\":"<<margin<<",\"exact_branch_margin\":"<<(I(margin)-I(rho)).l<<",\"charges\":[[0,0,1],[12,0,-1]],\"candidate_energy\":["<<candidateEnergy.l<<","<<candidateEnergy.h<<"],\"saddle_energy\":["<<(candidateEnergy-I(energyError)).l<<","<<(candidateEnergy+I(energyError)).h<<"],\"energy_correction_abs_upper\":"<<energyError<<",\"candidate_relative_to_reference\":["<<relativeEnergy.l<<","<<relativeEnergy.h<<"],\"saddle_relative_to_reference\":["<<(relativeEnergy-I(energyError)).l<<","<<(relativeEnergy+I(energyError)).h<<"],\"morse_index\":1,\"symmetry_restriction_on_perturbations\":false}\n";
 report_guards();return 0;
}catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
