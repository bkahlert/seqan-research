//autor:Jakob

#ifndef GINGER_TEST_GSEARCH_GSCOREPEPTIDEQGRAM_H_
#define GINGER_TEST_GSEARCH_GSCOREPEPTIDEQGRAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch/GScorePeptideQGram.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GScorePeptideQGram_CQ_AA)
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
	
	int res=GScorePeptideQGram(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
	SEQAN_ASSERT_EQ(score,300);
}

SEQAN_DEFINE_TEST(test_GScorePeptideQGram_LA_AA)
{
	GFastaRecord<Peptide> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="SSSSSHPAGCIIIIIIIIIII";
	GFastaRecord<Peptide> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="LLLLLLLHPAGC";
	Score<int> scoringScheme(1, -2, -1);
	string scoreMode="LOCAL_ALIGNMENT_MODE";
	int score;
	
	int res=GScorePeptideQGram(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
}

SEQAN_DEFINE_TEST(test_GScorePeptideQGram_GA_AA)
{
	GFastaRecord<Peptide> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="SSSSSHPAGCIIIIIIIIIII";
	GFastaRecord<Peptide> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="LLLLLLLHPAGC";
	Score<int> scoringScheme(1, -2, -1);
	string scoreMode="GLOBAL_ALIGNMENT_MODE";
	int score;
	
	int res=GScorePeptideQGram(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
}


SEQAN_DEFINE_TEST(test_GScorePeptideQGram_CQ_AA_empty)
{
	GFastaRecord<Peptide> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="";
	GFastaRecord<Peptide> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="";
	Score<int> scoringScheme(1, -2, -1);
	string scoreMode="COMMON_QGRAM_MODE";
	int score;
	
	int res=GScorePeptideQGram(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
	SEQAN_ASSERT_EQ(score,0);
}

SEQAN_DEFINE_TEST(test_GScorePeptideQGram_LA_AA_empty)
{
	GFastaRecord<Peptide> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="";
	GFastaRecord<Peptide> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="";
	Score<int> scoringScheme(1, -2, -1);
	string scoreMode="LOCAL_ALIGNMENT_MODE";
	int score;
	
	int res=GScorePeptideQGram(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
	SEQAN_ASSERT_EQ(score,0);
}

SEQAN_DEFINE_TEST(test_GScorePeptideQGram_GA_AA_empty)
{
	GFastaRecord<Peptide> queryRecord;
	queryRecord.id="testQueryId";
	queryRecord.seq="";
	GFastaRecord<Peptide> clusterRecord;
	clusterRecord.id="testClusterId";
	clusterRecord.seq="";
	Score<int> scoringScheme(1, -2, -1);
	string scoreMode="GLOBAL_ALIGNMENT_MODE";
	int score;
	
	int res=GScorePeptideQGram(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
	SEQAN_ASSERT_EQ(score,0);
}


#endif  // GINGER_TEST_GSEARCH_GSCOREPEPTIDEQGRAM_H_ 
