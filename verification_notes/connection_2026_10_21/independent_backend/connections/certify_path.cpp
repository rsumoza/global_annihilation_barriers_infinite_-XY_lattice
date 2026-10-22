#define main original_chain_main
#include "engine.cpp"
#undef main

// Streaming validation: memory is independent of the number of path nodes.
int main(int argc,char**argv){try{
 need(argc>=2&&argc<=4,"usage: certify_path INPUT_PREFIX [OUTPUT_TABLE_PREFIX] [--leg|--terminal-leg]");bool leg=argc==4&&std::string(argv[3])=="--leg",terminal=argc==4&&std::string(argv[3])=="--terminal-leg";need(argc<4||leg||terminal,"unknown path mode");std::string pre=argv[1],table=argc>=3?argv[2]:pre;
 need(LDBL_MANT_DIG>=64&&std::numeric_limits<double>::is_iec559&&sizeof(double)==8,"platform");need(std::fegetround()==FE_TONEAREST,"rounding");uint16_t endian=1;need(*(unsigned char*)&endian==1,"endian");
 std::ifstream cfg(pre+".cfg");int nodes=0,N=0,L=0;need((bool)(cfg>>nodes>>N>>L),"path config");need(nodes>=2&&nodes<50000&&N>=32&&N<=512&&N%2==0&&L==12,"path domain");
 auto h=readraw(pre+".reference_h.raw",2048ull*2048);int start=1024-N/2;R H=0,Hfar=0;
 for(int i=0;i<2048;i++)for(int j=0;j<2048;j++){
  R v=h[(size_t)i*2048+j];H=std::max(H,std::fabs(v));need((i&&j&&i<2047&&j<2047)||v==0,"reference boundary");
  if(!(i>=start+1&&i<start+N-1&&j>=start+1&&j<start+N-1))Hfar=std::max(Hfar,std::fabs(v));
 }
 need(H<.2L,"reference phase domain");I pi=PI();I exteriorAngle=I(12)/sq(I(N/2-8))+I(2)*I(Hfar);need(exteriorAngle.h<pi.l,"outside charges or branch uncontrolled");
 std::vector<CS>base((size_t)N*N);for(int i=0;i<N;i++)for(int j=0;j<N;j++)base[(size_t)i*N+j]=spin(i-(N-1)/2.0L+6,j-(N-1)/2.0L,12,h[(size_t)(start+i)*2048+start+j]);
 auto initial=charge_patch(base,N,6,pi);need(initial.q.size()==2,"reference charge count");bool plus=false,minus=false;for(auto v:initial.q){plus|=v==std::array<int,3>{0,0,1};minus|=v==std::array<int,3>{12,0,-1};}need(plus&&minus,"reference charges");
 std::array<std::vector<CS>,2>bond;bond[0].resize((size_t)(N-1)*N);bond[1].resize((size_t)N*(N-1));
 for(int ax=0;ax<2;ax++)for(int i=0;i<N-(ax==0);i++)for(int j=0;j<N-(ax==1);j++)bond[ax][(size_t)i*(N-(ax==1))+j]=product(base[(size_t)(i+(ax==0))*N+j+(ax==1)],base[(size_t)i*N+j]);
 std::ifstream file(pre+".z.raw",std::ios::binary|std::ios::ate);need((bool)file,"path payload open");uint64_t bytes=(uint64_t)nodes*N*N*sizeof(double);need(file.tellg()==(std::streamoff)bytes,"path payload length");file.seekg(0);
 std::ofstream ee(table+".exact.hex");ee<<std::hexfloat;
 std::ofstream nodeTable(table+".nodes.jsonl"),segmentTable(table+".segments.jsonl");need((bool)nodeTable&&(bool)segmentTable,"table output");nodeTable<<std::setprecision(24);segmentTable<<std::setprecision(24);
 std::vector<double>z((size_t)N*N),prev((size_t)N*N);R prevEnergy=0,maximum=0,maxCurve=0;int worst=-1;I endpointEnergy,firstEnergy;R maximumGradientNorm=0;
 for(int k=0;k<nodes;k++){
  file.read((char*)z.data(),(std::streamsize)z.size()*sizeof(double));need((bool)file,"path node read");
  for(int i=0;i<N;i++)for(int j=0;j<N;j++){R v=z[(size_t)i*N+j];need(std::isfinite(v),"nonfinite node");need((i&&j&&i<N-1&&j<N-1)||v==0,"path node boundary");if(k==0&&!leg&&!terminal)need(v==0,"path first node not zero");}
  I energy(0);R curve=0,norm=0;
  for(int ax=0;ax<2;ax++)for(int i=0;i<N-(ax==0);i++)for(int j=0;j<N-(ax==1);j++){
   size_t v=(size_t)i*N+j,t=(size_t)(i+(ax==0))*N+j+(ax==1);I d=z[t]==z[v]?I(0):I(z[t])-I(z[v]);CS tr=mag(d)<=.25L?small(d):trig(d),b=bond[ax][(size_t)i*(N-(ax==1))+j];
   energy=energy+b.c*(I(1)-tr.c)+b.s*tr.s;norm=normsq_up_add(norm,d);
   if(k){I jump=(I(z[t])-I(prev[t]))-(I(z[v])-I(prev[v]));curve=normsq_up_add(curve,jump);}
  }
  ee<<"NODE "<<k<<" "<<energy.l<<" "<<energy.h<<" "<<normup(norm)<<"\n";
 maximumGradientNorm=std::max(maximumGradientNorm,normup(norm));nodeTable<<"{\"node\":"<<k<<",\"relative_energy\":["<<energy.l<<","<<energy.h<<"],\"gradient_norm_upper\":"<<normup(norm)<<"}\n";
  if(k){R curv=(I(curve)/I(8)).h,upper=(I(std::max(energy.h,prevEnergy))+I(curv)).h;maxCurve=std::max(maxCurve,curv);if(upper>maximum){maximum=upper;worst=k-1;}ee<<"SEG "<<k-1<<" "<<curv<<" "<<upper<<"\n";segmentTable<<"{\"segment\":"<<k-1<<",\"curvature_error_upper\":"<<curv<<",\"relative_energy_upper\":"<<upper<<"}\n";}
  if(k==0)firstEnergy=energy;prevEnergy=energy.h;endpointEnergy=energy;prev.swap(z);
  if(k%200==0)std::cerr<<"Validated path node "<<k<<"/"<<nodes<<"\n";
 }
 if(leg){
  nodeTable.close();segmentTable.close();need((bool)nodeTable&&(bool)segmentTable,"table close");
  ee<<"U path_relative_reference_upper "<<maximum<<"\n";
 std::cout<<std::setprecision(24)<<"{\"status\":\"COMPACT_PATH_LEG_PASS_RELATIVE_REFERENCE\",\"nodes\":"<<nodes<<",\"side\":"<<N<<",\"segments\":"<<nodes-1<<",\"path_relative_reference_upper\":"<<maximum<<",\"maximizing_segment\":"<<worst<<",\"maximum_curvature_error\":"<<maxCurve<<",\"initial_relative_energy\":["<<firstEnergy.l<<","<<firstEnergy.h<<"],\"endpoint_relative_energy\":["<<endpointEnergy.l<<","<<endpointEnergy.h<<"],\"all_time_intervals_covered\":true,\"minimum_bridge_dependency_verified_here\":false,\"endpoint_annihilation_verified_here\":false,\"passes_through_exact_saddle_verified_here\":false}\n";ee.close();need((bool)ee,"path exact sidecar output");report_guards();return 0;
 }
 // prev is the last validated node. Its exact unit phases determine charges.
 std::vector<CS>end((size_t)N*N);std::vector<I>theta((size_t)N*N);
 for(int i=0;i<N;i++)for(int j=0;j<N;j++){
  size_t v=(size_t)i*N+j;I zz(prev[v]);CS r=mag(zz)<=.25L?small(zz):trig(zz),b=base[v];end[v]={b.c*r.c-b.s*r.s,b.s*r.c+b.c*r.s};R x=i-(N-1)/2.0L+6,y=j-(N-1)/2.0L;
  theta[v]=angle(I(x),I(y),pi)-angle(I(x-12),I(y),pi)+I(h[(size_t)(start+i)*2048+start+j])+zz;
 }
 auto charges=charge_patch(end,N,6,pi);need(charges.q.empty(),"endpoint still has charges");
 std::array<std::vector<int>,2>slip;slip[0].resize((size_t)(N-1)*N);slip[1].resize((size_t)N*(N-1));R maximumPrincipal=0;
 for(int ax=0;ax<2;ax++)for(int i=0;i<N-(ax==0);i++)for(int j=0;j<N-(ax==1);j++){
  size_t s=(size_t)i*N+j,t=(size_t)(i+(ax==0))*N+j+(ax==1);I beta=principal(product(end[t],end[s]),pi);maximumPrincipal=std::max(maximumPrincipal,mag(beta));I turn=(theta[t]-theta[s]-beta)/(I(2)*pi);
  int winding=(int)std::llround((turn.l+turn.h)/2);need(turn.l>winding-.25L&&turn.h<winding+.25L,"integer winding not isolated");slip[ax][(size_t)i*(N-(ax==1))+j]=winding;
 }
 R margin=std::min((pi-I(maximumPrincipal)).l,(pi-exteriorAngle).l);need(margin>0,"endpoint principal branch not strict");
 std::vector<int>lift((size_t)N*N,0);for(int i=0;i<N;i++)for(int j=0;j<N;j++){
  if(j)lift[(size_t)i*N+j]=lift[(size_t)i*N+j-1]+slip[1][(size_t)i*(N-1)+j-1];else if(i)lift[(size_t)i*N]=lift[(size_t)(i-1)*N]+slip[0][(size_t)(i-1)*N];
 }
 int maxTurn=0;for(int i=0;i<N;i++)for(int j=0;j<N;j++){
  int n=lift[(size_t)i*N+j];maxTurn=std::max(maxTurn,std::abs(n));if(i==0||j==0||i==N-1||j==N-1)need(n==0,"integer lift not compact");
  if(i+1<N)need(lift[(size_t)(i+1)*N+j]-n==slip[0][(size_t)i*N+j],"horizontal integer lift mismatch");if(j+1<N)need(lift[(size_t)i*N+j+1]-n==slip[1][(size_t)i*(N-1)+j],"vertical integer lift mismatch");
 }
 nodeTable.close();segmentTable.close();need((bool)nodeTable&&(bool)segmentTable,"table close");
 ee<<"U path_relative_reference_upper "<<maximum<<"\n";
 std::cout<<std::setprecision(24)<<"{\"status\":\""<<(terminal?"COMPACT_TERMINAL_PATH_PASS_RELATIVE_REFERENCE":"COMPACT_ANNIHILATION_PATH_PASS_RELATIVE_REFERENCE")<<"\",\"nodes\":"<<nodes<<",\"side\":"<<N<<",\"segments\":"<<nodes-1<<",\"path_relative_reference_upper\":"<<maximum<<",\"maximizing_segment\":"<<worst<<",\"maximum_curvature_error\":"<<maxCurve<<",\"maximum_gradient_norm\":"<<maximumGradientNorm<<",\"endpoint_relative_energy\":["<<endpointEnergy.l<<","<<endpointEnergy.h<<"],\"endpoint_charge_count\":0,\"endpoint_branch_margin\":"<<margin<<",\"compact_integer_lift_max\":"<<maxTurn<<",\"monotone_continuation_to_constant\":true,\"all_time_intervals_covered\":true,\"minimum_bridge_dependency_verified_here\":false,\"passes_through_exact_saddle_verified_here\":false}\n";
 ee.close();need((bool)ee,"path exact sidecar output");report_guards();return 0;
}catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
