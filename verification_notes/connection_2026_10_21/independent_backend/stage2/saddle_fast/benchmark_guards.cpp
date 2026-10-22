#define main unused_main
#include "engine.cpp"
#undef main
#include <chrono>
int main(){I x(1),y(0.00001L);auto t=std::chrono::steady_clock::now();for(int i=0;i<2000000;i++)x=x+y;std::cout<<std::setprecision(30)<<x.l<<" "<<x.h<<" "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<"\n";report_guards();}
