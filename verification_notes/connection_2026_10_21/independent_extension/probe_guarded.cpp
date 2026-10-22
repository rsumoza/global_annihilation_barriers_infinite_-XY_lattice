#define main ignored_original_engine_main
#include "instrumented/engine.cpp"
#undef main
int main(){
 bool nan_rejected=false,reversed_rejected=false;
 try{positive_by_congruence({I(std::numeric_limits<R>::quiet_NaN())},{1.0},1);}catch(const std::exception&){nan_rejected=true;}
 try{I reversed(1,0);}catch(const std::exception&){reversed_rejected=true;}
 need(nan_rejected&&reversed_rejected,"guard negative control accepted");
 std::cout<<"NAN_MATRIX_REJECTED "<<nan_rejected<<"\nREVERSED_INTERVAL_REJECTED "<<reversed_rejected<<"\n";
 return 0;
}
