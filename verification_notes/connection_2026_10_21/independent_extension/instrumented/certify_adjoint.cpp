#define main original_engine_main
#include "engine.cpp"
#undef main

static CS spin_full(R x,R y,R hh){
 I X(x),Y(y),Z(x-12),den=root(sq(X)+sq(Y))*root(sq(Z)+sq(Y));
 CS p{(X*Z+sq(Y))/den,I(-12)*Y/den},r=trig(I(hh));
 return {p.c*r.c-p.s*r.s,p.s*r.c+p.c*r.s};
}
static void ji(I a){exact_interval(a);std::cout<<"[\""<<dn(a.l)<<"\",\""<<up(a.h)<<"\"]";}
int main(int argc,char**argv){try{
 need(argc==6,"usage: certify_adjoint side Qold Qnew Hsaddle REFLECT_0_OR_1");
 need(LDBL_MANT_DIG>=64&&sizeof(double)==8&&std::numeric_limits<double>::is_iec559,"platform");
 need(std::fegetround()==FE_TONEAREST,"rounding");
 int n=std::stoi(argv[1]),ref=std::stoi(argv[5]);need(n==512||n==1024,"side");need(ref==0||ref==1,"reflection");int M=2*n,off=(n-64)/2,ho=(4096-n)/2;
 auto Q=readraw(argv[2],8192),q=readraw(argv[3],4096),h=readraw(argv[4],4096ull*4096);
 for(int i=0;i<64;i++)for(int j=0;j<64;j++)if(!i||!j||i==63||j==63){
  need(q[i*64+j]==0,"q compact boundary");
  need(Q[i*64+j]==0&&Q[4096+i*64+j]==0,"Q compact boundary");
 }
 auto z=readraw("z"+std::to_string(n)+".raw",(size_t)n*n),u=readraw("u"+std::to_string(M)+".raw",(size_t)M*M);
 if(ref){
  for(int i=0;i<n/2;i++)for(int j=0;j<n;j++)std::swap(z[(size_t)i*n+j],z[(size_t)(n-1-i)*n+j]);
  for(int i=0;i<M/2;i++)for(int j=0;j<M;j++)std::swap(u[(size_t)i*M+j],u[(size_t)(M-1-i)*M+j]);
  for(int i=0;i<32;i++)for(int j=0;j<64;j++)std::swap(q[i*64+j],q[(63-i)*64+j]);
 }
 std::array<I,4>G{};std::array<I,2>g{};
 for(int ax=0;ax<2;ax++)for(int i=0;i<64-(ax==0);i++)for(int j=0;j<64-(ax==1);j++){
  int a=i*64+j,b=(i+(ax==0))*64+j+(ax==1);I d0=I(Q[b])-I(Q[a]),d1=I(Q[4096+b])-I(Q[4096+a]),dq=I(q[b])-I(q[a]);
  G[0]=G[0]+sq(d0);G[1]=G[1]+d0*d1;G[3]=G[3]+sq(d1);g[0]=g[0]+d0*dq;g[1]=g[1]+d1*dq;
 }
 G[2]=G[1];I det=G[0]*G[3]-sq(G[1]);need(det.l>0,"Gram");
 std::array<I,4>Gi{G[3]/det,-G[1]/det,-G[1]/det,G[0]/det};
 std::array<I,2>c{Gi[0]*g[0]+Gi[1]*g[1],Gi[2]*g[0]+Gi[3]*g[1]};I prnorm=root(g[0]*c[0]+g[1]*c[1]);c[0]=-c[0]/prnorm;c[1]=-c[1]/prnorm;
 std::vector<I>r((size_t)n*n),f((size_t)M*M);std::vector<CS>sp((size_t)n*n);
 for(int i=0;i<64;i++)for(int j=0;j<64;j++)r[(i+off)*n+j+off]=c[0]*I(Q[i*64+j])+c[1]*I(Q[4096+i*64+j]);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){
  need((i&&j&&i<n-1&&j<n-1)||z[(size_t)i*n+j]==0,"compact z boundary");
  sp[(size_t)i*n+j]=spin_full(i-(n-1)/2.L+6,j-(n-1)/2.L,h[(size_t)(ref?4095-i-ho:i+ho)*4096+j+ho]);
 }
 I pairing(0),zn2(0);R dzmax=0,sindzmax=0;
 for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){
  int a=i*n+j,b=(i+(ax==0))*n+j+(ax==1);I dz=I(z[b])-I(z[a]),dr=r[b]-r[a];CS cc=product(sp[b],sp[a]);
  I flux=dr-cc.c*dz;size_t A=(size_t)(i+n/2)*M+j+n/2,B=(size_t)(i+(ax==0)+n/2)*M+j+(ax==1)+n/2;
  f[A]=f[A]-flux;f[B]=f[B]+flux;pairing=pairing+cc.s*dz;zn2=zn2+sq(dz);
  dzmax=std::max(dzmax,mag(dz));sindzmax=std::max(sindzmax,mag(cc.s*dz));
 }
 // Exact f is a compact divergence, hence has zero mean. A finite Neumann
 // Poisson energy is an upper bound for its full-lattice dual energy.
 R en=0,res=0;I fu(0),mean(0);
 for(int i=0;i<M;i++)for(int j=0;j<M;j++){
  size_t a=(size_t)i*M+j;
  if(i+1<M)en=normsq_up_add(en,I(u[a+M])-I(u[a]));
  if(j+1<M)en=normsq_up_add(en,I(u[a+1])-I(u[a]));
  res=normsq_up_add(res,f[a]-lap_at(u.data(),M,i,j));
  fu=fu+f[a]*I(u[a]);mean=mean+f[a];
 }
 need(mean.l<=0&&mean.h>=0,"force interval mean");
 I eta=I(normup(en))+I(M)/I(2)*I(normup(res));
 // |Dh_e| <= ||h||X/sqrt(2), the exact nearest-neighbor resistance.
 // Taylor's integral form retains sin(D shat) before absolute values.
 I rho(.012825L),nonlinear=(I(sindzmax)+rho/root(I(2))*I(dzmax))*sq(rho)/I(2);
 I radius=eta*rho+nonlinear,total=I(mag(pairing))+radius,correction=-pairing+I(-radius.h,radius.h);
 need(total.h<.001L,"offset threshold");
 exact_constant("rho_safe",.012825L);
 std::cout<<std::setprecision(30)<<"{\"status\":\"ADJOINT_NORMAL_OFFSET_PASS\",\"side\":"<<n<<",\"neumann_side\":"<<M<<",\"reflected_x\":"<<ref<<",\"rho_safe\":\".012825\",\"r_coefficients\":[";ji(c[0]);std::cout<<",";ji(c[1]);
 std::cout<<"],\"force_pairing\":";ji(pairing);std::cout<<",\"z_Xnorm\":";ji(root(zn2));
 std::cout<<",\"max_Dz_upper\":\""<<up(dzmax)<<"\",\"max_sin_Dshat_Dz_upper\":\""<<up(sindzmax)<<"\",\"neumann_poisson_residual_l2_upper\":\""<<normup(res)<<"\",\"adjoint_residual_Xdual_upper\":\""<<eta.h<<"\",\"nonlinear_remainder_upper\":\""<<nonlinear.h<<"\",\"normal_offset_error_upper\":\""<<total.h<<"\",\"signed_normal_correction\":";ji(correction);std::cout<<"}\n";
 finish_exact();return 0;
}catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
