#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include "PEMer_Lite.h"



SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_median)
{
	candidate test;
	modifier tester;
	char *file =" /Informatik/Development/test.sam";
	//saminput(test,file);
	SEQAN_ASSERT_EQ(50,50);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
   SEQAN_CALL_TEST(test_my_app_PEMer_Lite_median);
}
SEQAN_END_TESTSUITE