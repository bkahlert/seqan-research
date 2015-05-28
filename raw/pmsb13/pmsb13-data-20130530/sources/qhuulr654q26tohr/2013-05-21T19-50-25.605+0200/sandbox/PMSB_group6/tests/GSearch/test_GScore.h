#ifndef GINGER_TEST_GSEARCH_GSCORE_H_
#define GINGER_TEST_GSEARCH_GSCORE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch/GScore.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GScore_CQ_Dna)
{
	GFastaRecord<Dna5String> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="TTTATTGCTAAGCCAAA";
	GFastaRecord<Dna5String> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="CCCCATTGCTAAGCCTT";
	Score<int> scoringScheme(1, -2, -1);
	string scoreMode="COMMON_QGRAM_MODE";
	int score;
	
	int res=GScore(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
	SEQAN_ASSERT_EQ(score,400);
}

SEQAN_DEFINE_TEST(test_GScore_CQ_Rna)
{
	GFastaRecord<Rna5String> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="UUUAUUGCUAAGCCAAA";
	GFastaRecord<Rna5String> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="CCCCAUUGCUAAGCCUU";
	Score<int> scoringScheme(1, -2, -1);
	string scoreMode="COMMON_QGRAM_MODE";
	int score;
	
	int res=GScore(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
	SEQAN_ASSERT_EQ(score,400);
}

SEQAN_DEFINE_TEST(test_GScore_CQ_AA)
{
	GFastaRecord<Peptide> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="SSSSSHPAGCIIIIIIIIIII";
	GFastaRecord<Peptide> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="LLLLLLLHPAGC";
	Score<int> scoringScheme(1, -2, -1);
	string scoreMode="COMMON_QGRAM_MODE";
	int score;
	
	int res=GScore(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
	SEQAN_ASSERT_EQ(score,300);
}

#endif  // GINGER_TEST_GSEARCH_GSCORE_H_ 
