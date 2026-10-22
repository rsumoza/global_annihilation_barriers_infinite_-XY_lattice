#define main original_main
#include "baseline/engine.cpp"
#undef main
int main(){
 R nan=std::numeric_limits<R>::quiet_NaN();
 bool accepted=false;R v=0;try{v=positive_by_congruence({I(nan)},{1.0},1);accepted=true;}catch(const std::exception&e){std::cout<<"NAN_REJECT "<<e.what()<<"\n";}
 std::cout<<"NAN_ACCEPTED "<<accepted<<" returned "<<v<<"\n";
 R x=std::nextafter(1.0L,std::numeric_limits<R>::infinity());
 std::cout<<std::setprecision(22)<<"SERIAL_DEC "<<x<<"\n"<<std::hexfloat<<"SERIAL_HEX "<<x<<"\n"<<"CHART_DELTA_HEX "<<.001166548L<<"\n"<<"CHART_DELTA_PAD_UP_HEX "<<up(.001166548L)<<"\n"<<std::defaultfloat;
 return 0;
}
