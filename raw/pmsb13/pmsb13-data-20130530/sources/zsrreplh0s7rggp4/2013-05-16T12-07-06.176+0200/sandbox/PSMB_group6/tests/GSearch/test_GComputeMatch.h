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
  GFastaRecord<CharString> clusterRecord;
  TScoringScheme scoringScheme= ;
  TScore threshold= 300;
  string scoreMode="COMMON_QGRAM_MODE";

  
  int res=GComputeMatch(String<GMatch<TScore> > & matches, GFastaRecord<TSequence> & queryRecord, GFastaRecord<TSequence> & clusterRecord, TScoringScheme & scoringScheme, TScore threshold, string & scoreMode);
  SEQAN_ASSERT_EQ(res=0)
}

#endif  // GINGER_TEST_GCOMPUTEMATCH_H_ 
