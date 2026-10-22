#define main saved_engine_main
#include "../xy_optimality_extension/extension/engine.cpp"
#undef main
int main(){
 try {
  std::cout<<std::setprecision(30);
  const R M=std::numeric_limits<R>::max()/4;
  std::vector<I> A{I(M),I(M),I(M),I(-M)};
  std::vector<double> V{8,8,8,-8};
  std::cout<<"finite_input="<<std::isfinite(M)<<" indefinite_finite_input_test=";
  try {R p=positive_by_congruence(A,V,2);std::cout<<"ACCEPT margin="<<p<<"\n";}catch(const std::exception&e){std::cout<<"REJECT "<<e.what()<<"\n";}
  std::cout<<"explicit_nan_test=";
  try {R p=positive_by_congruence({I(std::numeric_limits<R>::quiet_NaN())},{1.0},1);std::cout<<"ACCEPT margin="<<p<<"\n";}catch(const std::exception&e){std::cout<<"REJECT "<<e.what()<<"\n";}
  auto x=trig(I(-7.99L,7.99L));
  std::cout<<"wide_trig_finite="<<(std::isfinite(x.c.l)&&std::isfinite(x.c.h)&&std::isfinite(x.s.l)&&std::isfinite(x.s.h))<<" cos="<<x.c.l<<","<<x.c.h<<" sin="<<x.s.l<<","<<x.s.h<<"\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}
}
