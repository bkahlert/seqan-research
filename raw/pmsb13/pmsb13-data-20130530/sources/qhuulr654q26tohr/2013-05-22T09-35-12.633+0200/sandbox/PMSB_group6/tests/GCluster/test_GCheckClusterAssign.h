//autor:Jakob

#ifndef GINGER_TEST_GCHECKCLUSTERASSIGN_H_
#define GINGER_TEST_GCHECKCLUSTERASSIGN_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GCluster/GCheckClusterAssign.h>

using namespace seqan;
using namespace std;


SEQAN_DEFINE_TEST(test_GCheckClusterAssign_dummy)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCAAA";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACC";
	
	int score = 10;
	int threshold = 5;
	double lengthThreshold=0.5;
	
	SEQAN_ASSERT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}

SEQAN_DEFINE_TEST(test_GCheckClusterAssign_lengthExact1)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCA";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACC";
	
	int score = 10;
	int threshold = 5;
	double lengthThreshold=0.5;
	
	SEQAN_ASSERT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}

SEQAN_DEFINE_TEST(test_GCheckClusterAssign_lengthExact2)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCAAAGTTTAGTAGT";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACC";
	
	int score = 10;
	int threshold = 5;
	double lengthThreshold=0.5;
	
	SEQAN_ASSERT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}

SEQAN_DEFINE_TEST(test_GCheckClusterAssign_lengthTooLong)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCAAAGTTTAGTAGTC";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACC";
	
	int score = 10;
	int threshold = 5;
	double lengthThreshold=0.5;
	
	SEQAN_ASSERT_NOT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}

SEQAN_DEFINE_TEST(test_GCheckClusterAssign_lengthTooShort)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCC";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACC";
	
	int score = 10;
	int threshold = 5;
	double lengthThreshold=0.5;
	
	SEQAN_ASSERT_NOT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}

SEQAN_DEFINE_TEST(test_GCheckClusterAssign_ScoreExact)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCAAA";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACC";
	
	int score = 5;
	int threshold = 5;
	double lengthThreshold=0.5;
	
	SEQAN_ASSERT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}

SEQAN_DEFINE_TEST(test_GCheckClusterAssign_ScoreTooLow)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCAAA";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACC";
	
	int score = 4;
	int threshold = 5;
	double lengthThreshold=0.5;
	
	SEQAN_ASSERT_NOT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}


#endif  // GINGER_TEST_GCHECKCLUSTERASSIGN_H_
