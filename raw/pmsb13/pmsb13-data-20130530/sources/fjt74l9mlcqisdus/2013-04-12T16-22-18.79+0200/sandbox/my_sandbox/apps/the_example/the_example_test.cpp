#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>


double quadrat(double var){
	return (var*var);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_3)
{
    double x = quadrat(3.0);
	SEQAN_ASSERT_IN_DELTA(x,9.0,0.1);
	//SEQAN_FAIL("Hello, tester!");
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat12)
{
    double x = quadrat(1.2);
	SEQAN_ASSERT_IN_DELTA(x,1.44,0.1);
	//SEQAN_FAIL("Hello, tester!");
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_funcs_quadrat3);
	SEQAN_CALL_TEST(test_my_app_funcs_quadrat12);
}
SEQAN_END_TESTSUITE