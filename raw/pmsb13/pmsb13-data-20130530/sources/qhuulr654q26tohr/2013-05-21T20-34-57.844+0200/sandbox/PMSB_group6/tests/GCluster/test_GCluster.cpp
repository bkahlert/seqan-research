#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_GCluster_systemTest.h"


SEQAN_BEGIN_TESTSUITE(test_GCluster)
{
	// Call tests.
	SEQAN_CALL_TEST(test_GCluster_systemTest_emptyFiles_LAMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_emptyFiles_CQMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_toyFiles_LAMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_toyFiles_CQMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_rfamFiles_LAMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_rfamFiles_GAMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_rfamFiles_CQMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_pfamFiles_LAMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_pfamFiles_GAMode);
	SEQAN_CALL_TEST(test_GCluster_systemTest_pfamFiles_CQMode);
}
SEQAN_END_TESTSUITE
