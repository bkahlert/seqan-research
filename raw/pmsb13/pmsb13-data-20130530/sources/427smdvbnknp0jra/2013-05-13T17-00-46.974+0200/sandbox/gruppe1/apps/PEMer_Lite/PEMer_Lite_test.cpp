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
	
	SEQAN_ASSERT_EQ(50,50);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
   SEQAN_CALL_TEST(test_my_app_PEMer_Lite_median);
}
SEQAN_END_TESTSUITE