#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>
#include "read_stats.h"


SEQAN_DEFINE_TEST(test_stats_inits_score_count_to_zero){
    ReadStats stats(9);
    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[8][8],0);
    
}

SEQAN_BEGIN_TESTSUITE(test_seqc)
{
    SEQAN_CALL_TEST(test_stats_inits_score_count_to_zero);

}

SEQAN_END_TESTSUITE
