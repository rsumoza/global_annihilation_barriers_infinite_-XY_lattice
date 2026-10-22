#define main ignored_old_main
#include "engine.cpp"
#undef main
static CS spin_general(R x,R y,R hh){I X(x),Y(y),Z(x-12),den=root(sq(X)+sq(Y))*root(sq(Z)+sq(Y));CS p={(X*Z+sq(Y))/den,I(-12)*Y/den},r=trig(I(hh));return {p.c*r.c-p.s*r.s,p.s*r.c+p.c*r.s};}
static void jout(I a){exact_interval(a);std::cout<<"[\""<<dn(a.l)<<"\",\""<<up(a.h)<<"\"]";}
int main(int argc,char**argv){try{
 need(argc==6,"usage certify_projection Qold Qnew Hsaddle Hbase REFLECT_0_OR_1");int ref=std::stoi(argv[5]);need(ref==0||ref==1,"reflection");
 need(LDBL_MANT_DIG>=64&&sizeof(double)==8&&std::numeric_limits<double>::is_iec559,"platform");need(std::fegetround()==FE_TONEAREST,"rounding");
 constexpr int n=256,off=96,ho=1920,bo=896;
 auto qold=readraw(argv[1],8192),qnew=readraw(argv[2],4096),h=readraw(argv[3],4096ull*4096),base=readraw(argv[4],2048ull*2048),u=readraw("eta_u256.raw",n*n),ur=readraw("eta_r_u256.raw",n*n);
 if(ref)for(int i=0;i<n/2;i++)for(int j=0;j<n;j++){std::swap(u[i*n+j],u[(n-1-i)*n+j]);std::swap(ur[i*n+j],ur[(n-1-i)*n+j]);}auto hs=[&](int i,int j)->R{return h[(size_t)(ref?4095-i:i)*4096+j];};
 std::array<I,4>G{};std::array<I,2>g{},hproj{};I nq(0),Hqq(0);
 std::vector<I>q(n*n),Q0(n*n),Q1(n*n),f(n*n);std::vector<CS>sp(n*n);
 for(int i=0;i<64;i++)for(int j=0;j<64;j++){
  need((i&&j&&i<63&&j<63)||qnew[(ref?63-i:i)*64+j]==0,"qnew boundary");
  for(int k=0;k<2;k++)need((i&&j&&i<63&&j<63)||qold[k*4096+i*64+j]==0,"qold boundary");
  q[(i+off)*n+j+off]=I(qnew[(ref?63-i:i)*64+j]);Q0[(i+off)*n+j+off]=I(qold[i*64+j]);Q1[(i+off)*n+j+off]=I(qold[4096+i*64+j]);
 }
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)sp[i*n+j]=spin_general(i-(n-1)/2.L+6,j-(n-1)/2.L,hs(i+ho,j+ho));
 for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){
  int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);I dq=q[t]-q[v],d0=Q0[t]-Q0[v],d1=Q1[t]-Q1[v];
  G[0]=G[0]+sq(d0);G[1]=G[1]+d0*d1;G[3]=G[3]+sq(d1);g[0]=g[0]+d0*dq;g[1]=g[1]+d1*dq;nq=nq+sq(dq);Hqq=Hqq+product(sp[t],sp[v]).c*sq(dq);
  I dz=(I(hs(i+ho+(ax==0),j+ho+(ax==1)))-I(hs(i+ho,j+ho)))-(I(base[(size_t)(i+bo+(ax==0))*2048+j+bo+(ax==1)])-I(base[(size_t)(i+bo)*2048+j+bo]));
  hproj[0]=hproj[0]+d0*dz;hproj[1]=hproj[1]+d1*dz;
 }
 G[2]=G[1];I det=G[0]*G[3]-sq(G[1]);need(det.l>0,"Gram");std::array<I,4>Gi{G[3]/det,-G[1]/det,-G[1]/det,G[0]/det};std::array<I,2>c{Gi[0]*g[0]+Gi[1]*g[1],Gi[2]*g[0]+Gi[3]*g[1]};I proj=g[0]*c[0]+g[1]*c[1],d2=I(1)-proj/nq;need(d2.l>0&&d2.h<.5L,"delta");I nu=Hqq/nq,norm=root(nq),prnorm=root(proj);
 for(auto&z:q)z=z/norm;
 for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){
  int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);I z=(product(sp[t],sp[v]).c-nu)*(q[t]-q[v]);f[v]=f[v]-z;f[t]=f[t]+z;
 }
 I mean(0),fu(0);R rsq=0,fsq=0;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){int v=i*n+j;mean=mean+f[v];fu=fu+f[v]*I(u[v]);rsq=normsq_up_add(rsq,f[v]-lap_at(u.data(),n,i,j));fsq=normsq_up_add(fsq,f[v]);}
 need(mean.l<=0&&mean.h>=0,"mean");R l2err=(I(128)*I(128)*I(normup(rsq))).h;I fnorm=I(normup(fsq)),eta2=fu+I(-1,1)*fnorm*I(l2err);need(eta2.l>0&&l2err<1e-7L,"Poisson residual");I eta=root(eta2);
 std::vector<I>rv(n*n),fr(n*n);for(int v=0;v<n*n;v++)rv[v]=(c[0]*Q0[v]+c[1]*Q1[v])/prnorm;
 I nur(0);for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);nur=nur+product(sp[t],sp[v]).c*sq(rv[t]-rv[v]);}
 for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);I z=(product(sp[t],sp[v]).c-nur)*(rv[t]-rv[v]);fr[v]=fr[v]-z;fr[t]=fr[t]+z;}
 I meanr(0),fur(0);R rsqr=0,fsqr=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++){int v=i*n+j;meanr=meanr+fr[v];fur=fur+fr[v]*I(ur[v]);rsqr=normsq_up_add(rsqr,fr[v]-lap_at(ur.data(),n,i,j));fsqr=normsq_up_add(fsqr,fr[v]);}
 need(meanr.l<=0&&meanr.h>=0,"rforce mean");R l2errr=(I(128)*I(128)*I(normup(rsqr))).h;I eta2r=fur+I(-1,1)*I(normup(fsqr))*I(l2errr);need(eta2r.l>0&&l2errr<1e-7L,"rPoisson residual");I etar=root(eta2r);

 // rho=0.012825 strictly exceeds exact saddle radius. Candidate margin .0719549999 minus rho/sqrt(2) exceeds .062886.
 I rho(.012825L),m(.062886L),pert=rho/root(I(2));I A=-nu+pert,et=eta+pert;
 I mu=m*(I(1)-d2)-I(2)*et*root(d2*(I(1)-d2))-A*d2;need(mu.l>0,"projected-normal compression");I mucandidate=I(.0719549999L)*(I(1)-d2)-I(2)*eta*root(d2*(I(1)-d2))+nu*d2;I mustrong=mucandidate-pert;need(mustrong.l>mu.l,"improved single-charge compression");
 std::array<I,2>normal{-g[0]/prnorm,-g[1]/prnorm},as{Gi[0]*hproj[0]+Gi[1]*hproj[1],Gi[2]*hproj[0]+Gi[3]*hproj[1]};I offset=normal[0]*as[0]+normal[1]*as[1];std::array<I,2>tv{-normal[1],normal[0]};I tn=root(G[0]*sq(tv[0])+I(2)*G[1]*tv[0]*tv[1]+G[3]*sq(tv[1]));std::array<I,2>tangent{(G[0]*tv[0]+G[1]*tv[1])/tn,(G[2]*tv[0]+G[3]*tv[1])/tn};I tangentoffset=tangent[0]*as[0]+tangent[1]*as[1];
 exact_constant("rho_safe",.012825L);exact_constant("m_safe",.062886L);exact_constant("candidate_complement",.0719549999L);
 std::cout<<std::setprecision(30)<<"{\"status\":\"PROJECTED_NORMAL_POSITIVE_CAP_PASS\",\"reflected_x\":"<<ref<<",\"rho_safe\":\".012825\",\"positive_q_complement_safe\":\".062886\",\"q_raw_norm_squared\":";jout(nq);std::cout<<",\"delta_squared\":";jout(d2);std::cout<<",\"delta\":";jout(root(d2));std::cout<<",\"nu\":";jout(nu);std::cout<<",\"eta_candidate_upper\":\""<<up(eta.h)<<"\",\"poisson_l2_error_upper\":\""<<up(l2err)<<"\",\"mu_projected_normal\":";jout(mu);std::cout<<",\"candidate_plane_coercivity\":";jout(mucandidate);std::cout<<",\"mu_single_perturbation\":";jout(mustrong);std::cout<<",\"strong_cap_radius_threshold\":";jout(I(3)*root(I(2))*mustrong);std::cout<<",\"r_nu_candidate\":";jout(nur);std::cout<<",\"r_eta_candidate_upper\":\""<<up(etar.h)<<"\",\"r_mixed_true_upper\":\""<<up((etar+pert).h)<<"\",\"r_negative_true_upper\":\""<<up((-nur+pert).h)<<"\",\"r_poisson_l2_error_upper\":\""<<up(l2errr)<<"\",\"normal_coordinates\":[";jout(normal[0]);std::cout<<",";jout(normal[1]);std::cout<<"],\"tangent_coordinates\":[";jout(tangent[0]);std::cout<<",";jout(tangent[1]);std::cout<<"],\"candidate_tangent_offset\":";jout(tangentoffset);std::cout<<",\"candidate_saddle_coordinates\":[";jout(as[0]);std::cout<<",";jout(as[1]);std::cout<<"],\"candidate_plane_offset\":";jout(offset);std::cout<<",\"exact_plane_offset\":";jout(offset+I(-rho.h,rho.h));std::cout<<",\"positive_cap_radius_upper_threshold\":";jout(I(3)*root(I(2))*mu);std::cout<<",\"R_01_rim_gap_lower\":\""<<dn((mu*sq(I(.1L))/I(2)-pwr(I(.1L),3)/(I(6)*root(I(2)))).l)<<"\"}\n";
 finish_exact();return 0;
}catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
