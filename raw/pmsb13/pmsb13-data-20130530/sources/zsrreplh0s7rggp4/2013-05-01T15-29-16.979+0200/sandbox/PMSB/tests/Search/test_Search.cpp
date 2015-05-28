#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_Search_systemTest.h"


SEQAN_BEGIN_TESTSUITE(test_Search)
{
	// Call tests.
	SEQAN_CALL_TEST(test_Search_systemTest_emptyFiles);
	SEQAN_CALL_TEST(test_Search_systemTest_toyFiles);
}
SEQAN_END_TESTSUITE
