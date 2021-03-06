#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_cluster_align.h"
#include "test_cluster_assignCluster.h"
#include "test_cluster_base.h"
#include "test_cluster_systemTest.h"

SEQAN_BEGIN_TESTSUITE(test_cluster)
{
	// Call tests.
	SEQAN_CALL_TEST(test_cluster_align);
	SEQAN_CALL_TEST(test_cluster_assignCluster);
	SEQAN_CALL_TEST(test_cluster_base);
	SEQAN_CALL_TEST(test_cluster_systemTest_emptyFiles);
	SEQAN_CALL_TEST(test_cluster_systemTest_toyFiles);
}
SEQAN_END_TESTSUITE
