#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include "own_functions.h"

double quadrat(double var){
	return (var*var);
}

void iota(seqan::String<int> & result,int begin,int end){
	resize(result,end-begin,0);
	for (int i=begin,k=0;i<end;++k,++i){
		result[k]=i;
	}
}



SEQAN_DEFINE_TEST(test_my_app_hash)
{
	int result = hash(1,0,3);
	SEQAN_ASSERT_IN_DELTA(result,31.5,31.5);
	//SEQAN_FAIL("Hello, tester!");
}

SEQAN_DEFINE_TEST(test_my_app_funcs_iota2)
{
	seqan::String<int> result;
	iota(result,3,6);
	SEQAN_ASSERT_EQ(length(result),3u);
	SEQAN_ASSERT_EQ(result[0],3);
	SEQAN_ASSERT_EQ(result[1],4);
	SEQAN_ASSERT_EQ(result[2],5);
	//SEQAN_FAIL("Hello, tester!");
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat12)
{
    double x = quadrat(1.2);
	SEQAN_ASSERT_IN_DELTA(x,1.44,0.1);
	//SEQAN_FAIL("Hello, tester!");
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat3)
{
    double x = quadrat(3.0);
	SEQAN_ASSERT_IN_DELTA(x,9.0,0.1);
	//SEQAN_FAIL("Hello, tester!");
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_hash);
	//SEQAN_CALL_TEST(test_my_app_funcs_quadrat12);
	//SEQAN_CALL_TEST(test_my_app_funcs_iota);
	//SEQAN_CALL_TEST(test_my_app_funcs_iota2);

}
SEQAN_END_TESTSUITE