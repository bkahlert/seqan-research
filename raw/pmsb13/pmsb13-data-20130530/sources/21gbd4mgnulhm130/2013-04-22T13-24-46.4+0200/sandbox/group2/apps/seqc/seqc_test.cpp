#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>
#include "seqc.h"

double quadrat(double x){
    return x*x;
}

void iota(seqan::String<unsigned int> & result, unsigned int begin, unsigned int end){
   resize(result, end-begin);
    for(unsigned int i=0; i<length(result);i++){
       result[i]=i;
    }
}

SEQAN_DEFINE_TEST(test_stats_obj_init_counters)
{

    double x = quadrat(1.2);
    SEQAN_ASSERT_IN_DELTA(x,1.44,0.01); 
}

SEQAN_BEGIN_TESTSUITE(test_seqc)
{

    SEQAN_CALL_TEST(test_my_app_funcs_quadrat_12);
}

SEQAN_END_TESTSUITE
