#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_GSearch_systemTest.h"
#include "test_GCheckClusterMaster.h"
#include "test_GComputeMatch.h"
#include "test_GParseMasterId.h"
#include "test_GEvaluateScore.h"
#include "test_GPostProcessMatches.h"
#include "test_GWriteMatches.h"
SEQAN_BEGIN_TESTSUITE(test_Search)
{
	// Call tests.
	SEQAN_CALL_TEST(test_GParseMasterID_dummy);
	SEQAN_CALL_TEST(test_GParseMasterID_noMaster);
	SEQAN_CALL_TEST(test_GWriteMatches_dummy);
	SEQAN_CALL_TEST(test_GEvaluateScore_dummy);
	SEQAN_CALL_TEST(test_GEvaluateScore_dummy2);
	SEQAN_CALL_TEST(test_GEvaluateScore_sameId);
	SEQAN_CALL_TEST(test_GEvaluateScore_sameScore);
	SEQAN_CALL_TEST(test_GEvaluateScore_noId);
	SEQAN_CALL_TEST(test_GEvaluateScore_noScore);
	SEQAN_CALL_TEST(test_GPostProcessMatches_dummy1);
	SEQAN_CALL_TEST(test_GPostProcessMatches_dummy2);
	SEQAN_CALL_TEST(test_GPostProcessMatches_short1);
	SEQAN_CALL_TEST(test_GPostProcessMatches_short2);
	SEQAN_CALL_TEST(test_GPostProcessMatches_alreadySorted);
	SEQAN_CALL_TEST(test_GPostProcessMatches_empty);
	SEQAN_CALL_TEST(test_GComputeMatch_dummy);
	SEQAN_CALL_TEST(test_GComputeMatch_empty);
	SEQAN_CALL_TEST(test_GCheckClusterMaster_dummy);
	SEQAN_CALL_TEST(test_GCheckClusterMaster_halfMaster);
	SEQAN_CALL_TEST(test_GCheckClusterMaster_withoutMaster);
	SEQAN_CALL_TEST(test_GCheckClusterMaster_empty);
	SEQAN_CALL_TEST(test_GSearch_systemTest_emptyFiles_LAMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_emptyFiles_CQMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_emptyFiles_CQMode_unclustered);
	SEQAN_CALL_TEST(test_GSearch_systemTest_toyFiles_LAMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_toyFiles_CQMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_toyFiles_CQMode_unclustered);
	SEQAN_CALL_TEST(test_GSearch_systemTest_rfamFiles_LAMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_rfamFiles_LAMode_unclustered);
	SEQAN_CALL_TEST(test_GSearch_systemTest_rfamFiles_CQMode);
	SEQAN_CALL_TEST(test_GSearch_systemTest_speicherzugriffsfehler_0517);
	
}
SEQAN_END_TESTSUITE
