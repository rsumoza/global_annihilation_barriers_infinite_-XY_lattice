#define main saved_engine_main
#include "../xy_optimality_extension/extension/engine.cpp"
#undef main
int main(){
 R x=up(1.0L);
 std::cout<<"mantissa="<<LDBL_MANT_DIG<<"\n";
 std::cout<<"exact_hex="<<std::hexfloat<<x<<"\n";
 std::cout<<"bare_precision22="<<std::defaultfloat<<std::setprecision(22)<<x<<"\n";
 std::cout<<"padded_upper_precision30="<<std::setprecision(30)<<up(x)<<"\n";
 std::cout<<"padded_lower_precision30="<<std::setprecision(30)<<dn(x)<<"\n";
}
