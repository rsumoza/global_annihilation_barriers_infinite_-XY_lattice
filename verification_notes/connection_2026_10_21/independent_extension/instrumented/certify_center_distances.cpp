#define main ignored_old_main
#include "engine.cpp"
#undef main
static void jout(I a){exact_interval(a);std::cout<<"[\""<<dn(a.l)<<"\",\""<<up(a.h)<<"\"]";}
int main(int argc,char**argv){try{
 need(argc==6,"usage Qold Hsaddle Hbase centersdirectory REFLECT_0_OR_1");int ref=std::stoi(argv[5]);need(ref==0||ref==1,"reflection");need(LDBL_MANT_DIG>=64&&sizeof(double)==8&&std::numeric_limits<double>::is_iec559&&std::fegetround()==FE_TONEAREST,"platform");uint16_t en=1;need(*(unsigned char*)&en==1,"endian");
 constexpr int n=128,off=32,ho=1984,bo=960;
 auto qo=readraw(argv[1],8192),hs=readraw(argv[2],4096ull*4096),hb=readraw(argv[3],2048ull*2048);std::array<I,4>G{};std::array<I,2>dsproj{};I dnorm(0);
 auto diff=[&](int i,int j)->I{I z(hs[(size_t)(ref?4095-i:i)*4096+j]);if(i>=1024&&i<3072&&j>=1024&&j<3072)z=z-I(hb[(size_t)(i-1024)*2048+j-1024]);return z;};
 for(int ax=0;ax<2;ax++)for(int i=0;i<4096-(ax==0);i++)for(int j=0;j<4096-(ax==1);j++)dnorm=dnorm+sq(diff(i+(ax==0),j+(ax==1))-diff(i,j));
 for(int i=0;i<4096;i++)need(hs[i]==0&&hs[(size_t)4095*4096+i]==0&&hs[(size_t)i*4096]==0&&hs[(size_t)i*4096+4095]==0,"saddle boundary");
 for(int i=0;i<2048;i++)need(hb[i]==0&&hb[(size_t)2047*2048+i]==0&&hb[(size_t)i*2048]==0&&hb[(size_t)i*2048+2047]==0,"base boundary");
 std::array<std::vector<I>,2>q;for(auto&z:q)z.resize(n*n);for(int k=0;k<2;k++)for(int i=0;i<64;i++)for(int j=0;j<64;j++)q[k][(i+off)*n+j+off]=I(qo[k*4096+i*64+j]);
 for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){
  int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);I d0=q[0][t]-q[0][v],d1=q[1][t]-q[1][v],dd=diff(i+ho+(ax==0),j+ho+(ax==1))-diff(i+ho,j+ho);G[0]=G[0]+sq(d0);G[1]=G[1]+d0*d1;G[3]=G[3]+sq(d1);dsproj[0]=dsproj[0]+dd*d0;dsproj[1]=dsproj[1]+dd*d1;
 }
 G[2]=G[1];I det=G[0]*G[3]-sq(G[1]);need(det.l>0,"Gram");std::array<I,4>Gi{G[3]/det,-G[1]/det,-G[1]/det,G[0]/det};
 auto projectnorm=[&](std::array<I,2> g){return sq(g[0])*Gi[0]+I(2)*g[0]*g[1]*Gi[1]+sq(g[1])*Gi[3];};
 std::cout<<std::setprecision(30)<<"{\"status\":\"SADDLE_CENTER_Y_DISTANCES_PASS\",\"reflected_x\":"<<ref<<",\"candidate_difference_norm_squared\":";jout(dnorm);std::cout<<",\"candidate_Y_norm_squared\":";jout(dnorm-projectnorm(dsproj));std::cout<<",\"centers\":[";
 for(int k=0;k<109;k++){
  std::ostringstream fn;fn<<argv[4]<<"/center_"<<std::setw(3)<<std::setfill('0')<<k<<".raw";auto z=readraw(fn.str(),n*n);I zn(0),cross(0);std::array<I,2>zp{};
  for(int i=0;i<n;i++)need(z[i]==0&&z[(n-1)*n+i]==0&&z[i*n]==0&&z[i*n+n-1]==0,"center boundary");
  for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){
   int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);I dz=I(z[t])-I(z[v]),dd=diff(i+ho+(ax==0),j+ho+(ax==1))-diff(i+ho,j+ho);zn=zn+sq(dz);cross=cross+dz*dd;for(int l=0;l<2;l++)zp[l]=zp[l]+dz*(q[l][t]-q[l][v]);
  }
  I y2=dnorm+zn-I(2)*cross-projectnorm({dsproj[0]-zp[0],dsproj[1]-zp[1]});need(y2.h>=0,"negative Y squared upper");I yn=root(I(std::max((R)0,y2.l),y2.h));if(k)std::cout<<",";std::cout<<"{\"id\":"<<k<<",\"Y_distance_squared\":";jout(y2);std::cout<<",\"Y_distance\":";jout(yn);std::cout<<"}";
 }
 std::cout<<"]}\n";finish_exact();return 0;
}catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
