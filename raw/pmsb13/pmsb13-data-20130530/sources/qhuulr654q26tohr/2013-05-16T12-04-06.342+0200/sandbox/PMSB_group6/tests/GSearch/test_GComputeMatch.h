#ifndef GINGER_TEST_GSEARCH_GCOMPUTEMATCH_H_
#define GINGER_TEST_GSEARCH_GCOMPUTEMATCH_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GComputeMatch_dummy)
{
  String<GMatch<int> > matches;
  GFastaRecord<CharString> queryRecord;
  queryRecord.id="testQueryId";
  queryRecord.seq="DDDDtesttestQuerySeqDDDDtesttestQuerySeq";
  GFastaRecord<CharString> clusterRecord;
  clusterRecord.id="testClusterId";
  clusterRecord.seq="DDDDtesttestClusterSeqDDDDtesttestClusterSeq";
  Score<int> scoringScheme(1, -2, -1);
  int threshold= 1;
  string scoreMode="COMMON_QGRAM_MODE";

  int res=GComputeMatch(matches, queryRecord, clusterRecord, scoringScheme, threshold, scoreMode);
  SEQAN_ASSERT_EQ(res,0);
  SEQAN_ASSERT_EQ(length(matches),1);
  SEQAN_ASSERT_EQ(matches[0].id,"testClusterId");
  SEQAN_ASSERT_EQ(matches[0].queryId,"testQueryId");
  SEQAN_ASSERT_EQ(matches[0].score,400);
}

#endif  // GINGER_TEST_GCOMPUTEMATCH_H_ 
