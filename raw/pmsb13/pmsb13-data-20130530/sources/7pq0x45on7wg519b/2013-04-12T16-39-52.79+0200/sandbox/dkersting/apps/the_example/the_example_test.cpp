#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>

double quadrat(double x)
{
 	return x*x;
}

void iota(seqan::String<int> &result, int begin, int end)
{
	resize (result, end-begin);
	for(int i=0;i<length(result);++i)
		result[i]=i;
}

SEQAN_DEFINE_TEST(test_my_app_funcs_iota)
{
	seqan::String<int> result;
	iota(result,0,3);
	SEQAN_ASSERT_EQ(length(result),3);
	SEQAN_ASSERT_EQ(result[0],0);
	SEQAN_ASSERT_EQ(result[1],1);
	SEQAN_ASSERT_EQ(result[2],2);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_3)
{
    	double x = quadrat(3.0);
	SEQAN_ASSERT_EQ(x,9.0);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_12)
{
    	double x = quadrat(1.2);
	SEQAN_ASSERT_EQ(x,1.44);
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    	SEQAN_CALL_TEST(test_my_app_funcs_quadrat_3);
    	SEQAN_CALL_TEST(test_my_app_funcs_quadrat_12);
    	SEQAN_CALL_TEST(test_my_app_funcs_iota);
}
SEQAN_END_TESTSUITE
