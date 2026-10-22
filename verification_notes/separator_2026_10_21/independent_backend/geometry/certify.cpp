#define main ignored_old_main
#include "engine.cpp"
#undef main
struct CenterData {int id,gi,gj;std::string name;std::vector<double>z;std::array<I,2>g,a;I ynorm2;};
static void interval_out(I x){std::cout<<"[\""<<dn(x.l)<<"\",\""<<up(x.h)<<"\"]";}
int main(int argc,char**argv){try{
 need(argc==4,"usage: atlas_geometry Q64 CFG CENTER_DIRECTORY");
 need(LDBL_MANT_DIG>=64&&sizeof(double)==8&&std::numeric_limits<double>::is_iec559,"platform");need(std::fegetround()==FE_TONEAREST,"rounding");uint16_t endian=1;need(*(unsigned char*)&endian==1,"endian");
 auto rawq=readraw(argv[1],8192);std::ifstream cfg(argv[2]);int n=0,count=0;cfg>>n>>count;need((bool)cfg&&n==128&&count>=1&&count<1000,"geometry cfg");
 int off=(n-64)/2;std::array<std::vector<double>,2>q;for(auto&v:q)v.assign(n*n,0);
 for(int k=0;k<2;k++)for(int i=0;i<64;i++)for(int j=0;j<64;j++){need((i&&j&&i<63&&j<63)||rawq[k*4096+i*64+j]==0,"q boundary");q[k][(i+off)*n+j+off]=rawq[k*4096+i*64+j];}
 std::array<I,4>G{};for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){
  int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);I d0=I(q[0][t])-I(q[0][v]),d1=I(q[1][t])-I(q[1][v]);G[0]=G[0]+sq(d0);G[1]=G[1]+d0*d1;G[3]=G[3]+sq(d1);
 }G[2]=G[1];I det=G[0]*G[3]-sq(G[1]);need(det.l>0,"Gram inverse");std::array<I,4>Gi{G[3]/det,-G[1]/det,-G[1]/det,G[0]/det};
 std::vector<CenterData>C;for(int k=0;k<count;k++){
  CenterData c;cfg>>c.id>>c.gi>>c.gj>>c.name;need((bool)cfg&&c.id==k&&c.name.find('/')==std::string::npos&&c.name.find("..") ==std::string::npos,"center record");c.z=readraw(std::string(argv[3])+"/"+c.name,(size_t)n*n);I norm2(0);
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)need((i&&j&&i+1<n&&j+1<n)||c.z[i*n+j]==0,"center boundary");
  for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){
   int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);I dz=I(c.z[t])-I(c.z[v]);norm2=norm2+sq(dz);for(int l=0;l<2;l++)c.g[l]=c.g[l]+(I(q[l][t])-I(q[l][v]))*dz;
  }
  c.a={Gi[0]*c.g[0]+Gi[1]*c.g[1],Gi[2]*c.g[0]+Gi[3]*c.g[1]};c.ynorm2=norm2-c.g[0]*c.a[0]-c.g[1]*c.a[1];need(c.ynorm2.h>=0,"Y norm upper");C.push_back(std::move(c));
 }
 R maxDiff=0,maxY=0,maxCoordinateError=0;int pairs=0;std::vector<std::tuple<int,int,I>>distances;
 for(auto&c:C){maxY=std::max(maxY,root(I(0,c.ynorm2.h)).h);maxCoordinateError=std::max(maxCoordinateError,mag(c.a[0]-I(3)*I(c.gi)/I(25)));maxCoordinateError=std::max(maxCoordinateError,mag(c.a[1]-I(3)*I(c.gj)/I(25)));}
 for(int k=0;k<count;k++)for(int l=k+1;l<count;l++)if(std::abs(C[k].gi-C[l].gi)<=1&&std::abs(C[k].gj-C[l].gj)<=1){
  I norm2(0);for(int ax=0;ax<2;ax++)for(int i=0;i<n-(ax==0);i++)for(int j=0;j<n-(ax==1);j++){
   int v=i*n+j,t=(i+(ax==0))*n+j+(ax==1);I d=(I(C[k].z[t])-I(C[k].z[v]))-(I(C[l].z[t])-I(C[l].z[v]));norm2=norm2+sq(d);
  }
  std::array<I,2>g{C[k].g[0]-C[l].g[0],C[k].g[1]-C[l].g[1]};I yp=norm2-sq(g[0])*Gi[0]-I(2)*g[0]*g[1]*Gi[1]-sq(g[1])*Gi[3];need(yp.h>=0,"neighbor Y norm upper");I d=root(I(std::max((R)0,yp.l),yp.h));maxDiff=std::max(maxDiff,d.h);distances.emplace_back(k,l,d);pairs++;
 }
 std::ofstream ee("geometry_exact.hex");ee<<std::hexfloat;
 for(int k=0;k<4;k++)ee<<"G "<<k<<" "<<G[k].l<<" "<<G[k].h<<"\n";
 for(int k=0;k<4;k++)ee<<"Gi "<<k<<" "<<Gi[k].l<<" "<<Gi[k].h<<"\n";
 ee<<"MAX neighbor "<<maxDiff<<"\nMAX Y "<<maxY<<"\nMAX coord "<<maxCoordinateError<<"\n";
 for(int k=0;k<count;k++){auto&c=C[k];for(int a=0;a<2;a++)ee<<"A "<<k<<" "<<a<<" "<<c.a[a].l<<" "<<c.a[a].h<<"\n";ee<<"Y2 "<<k<<" "<<c.ynorm2.l<<" "<<c.ynorm2.h<<"\n";}
 for(auto[i,j,d]:distances)ee<<"D "<<i<<" "<<j<<" "<<d.l<<" "<<d.h<<"\n";
 ee.close();need((bool)ee,"exact geometry sidecar output");
 std::cout<<std::setprecision(30)<<"{\"status\":\"ATLAS_PROJECTION_AND_NEIGHBOR_GEOMETRY_PASS\",\"n\":"<<n<<",\"count\":"<<count<<",\"G\":[";
 for(int k=0;k<4;k++){if(k)std::cout<<",";interval_out(G[k]);}std::cout<<"],\"Gi\":[";for(int k=0;k<4;k++){if(k)std::cout<<",";interval_out(Gi[k]);}
 std::cout<<"],\"maximum_neighbor_Y_distance_upper\":\""<<up(maxDiff)<<"\",\"maximum_center_Y_norm_upper\":\""<<up(maxY)<<"\",\"maximum_nominal_coordinate_error_upper\":\""<<up(maxCoordinateError)<<"\",\"neighbor_pairs\":"<<pairs<<",\"centers\":[";
 for(int k=0;k<count;k++){if(k)std::cout<<",";auto&c=C[k];std::cout<<"{\"id\":"<<c.id<<",\"grid\":["<<c.gi<<","<<c.gj<<"],\"file\":\""<<c.name<<"\",\"a\":[";interval_out(c.a[0]);std::cout<<",";interval_out(c.a[1]);std::cout<<"],\"Y_norm_squared\":";interval_out(c.ynorm2);std::cout<<"}";}
 std::cout<<"],\"neighbors\":[";for(size_t k=0;k<distances.size();k++){if(k)std::cout<<",";auto[i,j,d]=distances[k];std::cout<<"{\"i\":"<<i<<",\"j\":"<<j<<",\"distance\":";interval_out(d);std::cout<<"}";}std::cout<<"]}\n";report_guards();return 0;
}catch(const std::exception&e){std::cerr<<"REJECT: "<<e.what()<<"\n";return 2;}}
