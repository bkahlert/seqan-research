#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>


double quadrat(double var){
	return (var*var);
}

void iota(seqan::String<int> & result,int begin,int end){
	resize(result,end-begin,0);
	for (int i=begin,k=0;i<end;++k,++i){
		result[k]=i;
	}
}



SEQAN_DEFINE_TEST(test_my_app_funcs_iota)
{
	seqan::String<int> result;
	iota(result,0,3);
	SEQAN_ASSERT_EQ(length(result),3);
	SEQAN_ASSERT_EQ(result[0],0);
	SEQAN_ASSERT_EQ(result[1],1);
	SEQAN_ASSERT_EQ(result[2],2);
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
	SEQAN_CALL_TEST(test_my_app_funcs_iota);
}
SEQAN_END_TESTSUITE