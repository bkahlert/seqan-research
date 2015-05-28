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
	clusterRecord.seq="TTGGGGGGGA";
	
	int score = 10;
	threshold = 5;
	lengthThreshold=0.5;
	
	SEQAN_ASSERT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}

SEQAN_DEFINE_TEST(test_GCheckClusterAssign_lengthExact)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCAAA";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACC";
	
	int score = 10;
	threshold = 5;
	lengthThreshold=0.5;
	
	SEQAN_ASSERT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}

SEQAN_DEFINE_TEST(test_GCheckClusterAssign_lengthTooLong)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCAAA";
	GFastaRecord<DnaString> clusterRecord;
	clusterRecord.id="myClusterMaster";
	clusterRecord.seq="TTGGGGGGGACCT";
	
	int score = 10;
	threshold = 5;
	lengthThreshold=0.5;
	
	SEQAN_ASSERT_NOT(GCheckClusterAssign(score, threshold, newRecord, clusterRecord, lengthThreshold));
}


#endif  // GINGER_TEST_GCHECKCLUSTERASSIGN_H_
