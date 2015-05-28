#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_clusterSearch.h"


SEQAN_BEGIN_TESTSUITE(test_clusterSearch)
{
	SEQAN_CALL_TEST(test_clusterSearch_emptyFiles);
}
SEQAN_END_TESTSUITE
