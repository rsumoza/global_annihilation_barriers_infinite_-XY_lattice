#define main original_chain_main
#include "engine.cpp"
#undef main
int main(int argc,char**argv){try{
 need(argc==7,"usage: certify_join_points FINE_H REFERENCE_H Q64 Z_MINUS Z_PLUS SIDE");int N=std::stoi(argv[6]);need(N>=32&&N<=512&&N%2==0,"point side");
 need(LDBL_MANT_DIG>=64&&std::numeric_limits<double>::is_iec559&&sizeof(double)==8,"platform");need(std::fegetround()==FE_TONEAREST,"rounding");uint16_t endian=1;need(*(unsigned char*)&endian==1,"endian");
 auto hf=readraw(argv[1],4096ull*4096),hb=readraw(argv[2],2048ull*2048),q=readraw(argv[3],4096);std::array<std::vector<double>,2>z{readraw(argv[4],(size_t)N*N),readraw(argv[5],(size_t)N*N)};
 for(int k=0;k<2;k++)for(int i=0;i<N;i++)for(int j=0;j<N;j++)need((i&&j&&i<N-1&&j<N-1)||z[k][(size_t)i*N+j]==0,"point support");
 for(int i=0;i<64;i++)for(int j=0;j<64;j++)need((i&&j&&i<63&&j<63)||q[i*64+j]==0,"q support");
 std::array<std::vector<I>,2>dq;dq[0].resize(63*64);dq[1].resize(64*63);I G(0);for(int ax=0;ax<2;ax++)for(int i=0;i<64-(ax==0);i++)for(int j=0;j<64-(ax==1);j++){
  I d=I(q[(i+(ax==0))*64+j+(ax==1)])-I(q[i*64+j]);dq[ax][i*(64-(ax==1))+j]=d;G=G+sq(d);
 }need(G.l>.5L&&G.h<2,"q scale");
 int fineStart=2048-N/2,baseStart=1024-N/2;
 auto difference=[&](int k,int i,int j){
  R f=(i>=0&&j>=0&&i<4096&&j<4096)?hf[(size_t)i*4096+j]:0;int bi=i-1024,bj=j-1024,zi=i-fineStart,zj=j-fineStart;
  R b=(bi>=0&&bj>=0&&bi<2048&&bj<2048)?hb[(size_t)bi*2048+bj]:0,w=(zi>=0&&zj>=0&&zi<N&&zj<N)?z[k][(size_t)zi*N+zj]:0;
  return I(b)+I(w)-I(f);
 };
 std::array<I,2>normsq{I(0),I(0)},projection{I(0),I(0)},energy{I(0),I(0)};
 for(int i=-1;i<4097;i++)for(int j=-1;j<4097;j++)for(int k=0;k<2;k++){
  I v=difference(k,i,j);if(i<4096)normsq[k]=normsq[k]+sq(difference(k,i+1,j)-v);if(j<4096)normsq[k]=normsq[k]+sq(difference(k,i,j+1)-v);
 }
 for(int ax=0;ax<2;ax++)for(int i=0;i<64-(ax==0);i++)for(int j=0;j<64-(ax==1);j++)for(int k=0;k<2;k++){
  int fi=2016+i,fj=2016+j;I d=difference(k,fi+(ax==0),fj+(ax==1))-difference(k,fi,fj);projection[k]=projection[k]+d*dq[ax][i*(64-(ax==1))+j];
 }
 std::vector<CS>base((size_t)N*N);for(int i=0;i<N;i++)for(int j=0;j<N;j++)base[(size_t)i*N+j]=spin(i-(N-1)/2.0L+6,j-(N-1)/2.0L,12,hb[(size_t)(baseStart+i)*2048+baseStart+j]);
 for(int ax=0;ax<2;ax++)for(int i=0;i<N-(ax==0);i++)for(int j=0;j<N-(ax==1);j++){
  size_t v=(size_t)i*N+j,t=(size_t)(i+(ax==0))*N+j+(ax==1);CS b=product(base[t],base[v]);for(int k=0;k<2;k++){I d=I(z[k][t])-I(z[k][v]);CS r=mag(d)<=.25L?small(d):trig(d);energy[k]=energy[k]+b.c*(I(1)-r.c)+b.s*r.s;}
 }
 std::ofstream ee("results/join_exact.hex");ee<<std::hexfloat<<"I q_gram "<<G.l<<" "<<G.h<<"\n";
 std::cout<<std::setprecision(24)<<"{\"status\":\"SADDLE_JOIN_POINT_COMPARISONS_PASS\",\"side\":"<<N<<",\"q_gram\":["<<G.l<<","<<G.h<<"],\"points\":[";
 for(int k=0;k<2;k++){
  I distance=root(I(std::max((R)0,normsq[k].l),normsq[k].h)),p=projection[k]/root(G),transverse=normsq[k]-sq(p);need(transverse.h>=0,"negative transverse square upper");R transverseUpper=root(I(0,std::max((R)0,transverse.h))).h;
  ee<<"I points/"<<k<<"/distance_to_saddle_candidate "<<distance.l<<" "<<distance.h<<"\nI points/"<<k<<"/signed_projection_on_normalized_q "<<p.l<<" "<<p.h<<"\nU points/"<<k<<"/transverse_norm_upper "<<transverseUpper<<"\nI points/"<<k<<"/relative_energy_to_reference "<<energy[k].l<<" "<<energy[k].h<<"\n";
 if(k)std::cout<<",";std::cout<<"{\"label\":\""<<(k?"plus":"minus")<<"\",\"distance_to_saddle_candidate\":["<<distance.l<<","<<distance.h<<"],\"signed_projection_on_normalized_q\":["<<p.l<<","<<p.h<<"],\"transverse_norm_upper\":"<<transverseUpper<<",\"relative_energy_to_reference\":["<<energy[k].l<<","<<energy[k].h<<"]}";
 }
 std::cout<<"],\"exact_saddle_error_not_applied_here\":true,\"hermite_join_verified_here\":false}\n";ee.close();need((bool)ee,"join exact sidecar output");report_guards();return 0;
}catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
