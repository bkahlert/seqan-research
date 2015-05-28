#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_structs.h"


SEQAN_BEGIN_TESTSUITE(test_structs)
{
	// Call tests.
	SEQAN_CALL_TEST(test_structs_FastaRecord);
	SEQAN_CALL_TEST(test_structs_ScoredSequence);
	SEQAN_CALL_TEST(test_structs_PerformanceSample);
}
SEQAN_END_TESTSUITE
