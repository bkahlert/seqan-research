#ifndef GINGER_TEST_GSEARCH_GCOMPUTEMATCH_H_
#define GINGER_TEST_GSEARCH_GCOMPUTEMATCH_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch/GComputeMatch.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GComputeMatch_dummy)
{
	String<GMatch<int, DnaString> > matches;
	GFastaRecord<DnaString> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="TTTATTGCTAAGCCAAA";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="CCCCATTGCTAAGCCTT";
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

SEQAN_DEFINE_TEST(test_GComputeMatch_empty)
{
	String<GMatch<int, DnaString> > matches;
	GFastaRecord<DnaString> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="TTTATTGCTAAGCCAAA";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="";
	clusterRecord.seq="";
	Score<int> scoringScheme(1, -2, -1);
	int threshold= 0;
	string scoreMode="COMMON_QGRAM_MODE";
	
	int res=GComputeMatch(matches, queryRecord, clusterRecord, scoringScheme, threshold, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
	SEQAN_ASSERT_EQ(length(matches),1);
	SEQAN_ASSERT_EQ(matches[0].id,"");
	SEQAN_ASSERT_EQ(matches[0].queryId,"testQueryId");
	SEQAN_ASSERT_EQ(matches[0].score,0);
}

#endif  // GINGER_TEST_GCOMPUTEMATCH_H_ 
