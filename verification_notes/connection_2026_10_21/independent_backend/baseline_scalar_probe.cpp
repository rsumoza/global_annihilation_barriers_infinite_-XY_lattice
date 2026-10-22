#define main original_main
#include "stage2/baseline/engine.cpp"
#undef main
int main(int argc,char**argv){try{
 need(argc==3,"usage: baseline_scalar_probe C_HEX RHO_HEX");R c=std::stold(argv[1]),rho=std::stold(argv[2]);R mu=(I(c)-I(rho)).l;
 std::cout<<std::hexfloat<<"RAW probe_c "<<c<<"\nRAW probe_rho "<<rho<<"\nL mu "<<mu<<"\nL branch_margin "<<1.12L<<"\nRAW compiled_rho "<<.017422L<<"\nRAW compiled_delta "<<.001166548L<<"\nRAW compiled_mustar "<<.05825L<<"\nRAW barrier_branch_limit "<<1.12L<<"\n";report_guards();return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}}
