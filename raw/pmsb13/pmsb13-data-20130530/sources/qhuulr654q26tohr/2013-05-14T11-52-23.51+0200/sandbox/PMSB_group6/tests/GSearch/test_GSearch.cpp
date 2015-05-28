#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_GSearch_systemTest.h"


SEQAN_BEGIN_TESTSUITE(test_Search)
{
	// Call tests.
	SEQAN_CALL_TEST(test_GSearch_systemTest_emptyFiles_LAMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_emptyFiles_CQMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_emptyFiles_CQMode_unclustered);
	SEQAN_CALL_TEST(test_GSearch_systemTest_toyFiles_LAMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_toyFiles_CQMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_toyFiles_CQMode_unclustered);
	SEQAN_CALL_TEST(test_GSearch_systemTest_rfamFiles_LAMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_rfamFiles_LAMode_unclustered);
	SEQAN_CALL_TEST(test_GSearch_systemTest_rfamFiles_CQMode);
}
SEQAN_END_TESTSUITE
