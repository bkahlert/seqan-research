#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include "PEMer_Lite.h"

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input)
{
	candidates test;
	modifier opt;
	statisticals stats
	char *file ="/Informatik/Development/test.sam";
	SEQAN_ASSERT(!saminput(test,opt,stats));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input2)
{
	candidate test;
	char *file =" ";
	SEQAN_ASSERT(saminput(test,file));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input3)
{
	candidate test;
	char *file ="2579816ß15362466183b513ü6ß148z5!§%%!&°%(!%§(%§c465cn4ü457346511425129562e85z148z5!§%%!&°%(!%§(%§&&&&&(";
	SEQAN_ASSERT(saminput(test,file));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_devide)
{
	candidate test,result;
	modifier opt;
	for(unsigned i=1;i<11;i++)
	{
		test.len.push_back(i);
		test.ref.push_back(i);
		test.pos.push_back(i);
		test.type.push_back(candidate::u);
	}
	opt.e=5;
	opt.sd=2;
	devide(test,result,opt);
	
	SEQAN_ASSERT_EQ(result.type[0],1);
	SEQAN_ASSERT_EQ(result.ref[0],1);
	SEQAN_ASSERT_EQ(result.type[1],1);
	SEQAN_ASSERT_EQ(result.ref[1],2);
	SEQAN_ASSERT_EQ(result.type[2],2);
	SEQAN_ASSERT_EQ(result.ref[2],8);
	SEQAN_ASSERT_EQ(result.type[3],2);
	SEQAN_ASSERT_EQ(result.ref[3],9);
	SEQAN_ASSERT_EQ(result.type[4],2);
	SEQAN_ASSERT_EQ(result.ref[4],10);
	
}


SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_median)
{
	candidate test;
	for(unsigned i=0;i<100;i++)
	{
		test.len.push_back(i);
	}
	median(test);
	SEQAN_ASSERT_EQ(test.e,50);
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_analyse)
{
	candidate test,result;

	for(unsigned i=0;i<100;i++)
	{
		test.len.push_back(i);
	}

	analyse(test);
	
	SEQAN_ASSERT_EQ(test.e,50);
	SEQAN_ASSERT_EQ(test.sd,28);
	
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_median);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input2);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input3);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_analyse);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_devide);

}


SEQAN_END_TESTSUITE