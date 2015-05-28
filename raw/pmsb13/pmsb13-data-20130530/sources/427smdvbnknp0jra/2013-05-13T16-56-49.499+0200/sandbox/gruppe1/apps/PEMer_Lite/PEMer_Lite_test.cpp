#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include "PEMer_Lite.h"



SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_median)
{
	candidate test;
	for(unsigned i=1;i<100;i++)
	{
		test.len.push_back(i);
	}
	
	SEQAN_ASSERT_EQ(median(test),500);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_3)
{
	double x =quadrat(3.0);
    SEQAN_ASSERT_IN_DELTA(x,9.0,0.01);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_12)
{
	double x =quadrat(1.2);
    SEQAN_ASSERT_IN_DELTA(x,1.44,0.01);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_funcs_quadrat_3);
	SEQAN_CALL_TEST(test_my_app_funcs_quadrat_12);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_median);
	SEQAN_CALL_TEST(test_my_app_funcs_iota2);
}
SEQAN_END_TESTSUITE