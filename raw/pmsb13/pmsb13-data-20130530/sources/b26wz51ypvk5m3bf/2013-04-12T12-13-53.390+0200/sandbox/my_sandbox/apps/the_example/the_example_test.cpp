#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/sequence.h>
using namespace seqan;
double quadrat(double d)
{
   return d*d;
}

SEQAN_DEFINE_TEST(test_my_app_funcs_iota)
{
    String<int> result;
    iota(result, 0, 3);
    SEQAN_ASSERT_EQ(length(result),3);
    SEQAN_ASSERT_EQ(result[0],0);
    SEQAN_ASSERT_EQ(result[1],1);
    SEQAN_ASSERT_EQ(result[2],2);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_von_drei)
{
   double x = quadrat(3.0); 
   SEQAN_ASSERT_IN_DELTA(x, 9.0,0.1);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_von_eins_komma_zwei)
{
   double x = quadrat(1.2); 
   SEQAN_ASSERT_IN_DELTA(x, 1.44,0.1);
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_funcs_quadrat_von_drei);
    SEQAN_CALL_TEST(test_my_app_funcs_quadrat_von_eins_komma_zwei);
}
SEQAN_END_TESTSUITE
