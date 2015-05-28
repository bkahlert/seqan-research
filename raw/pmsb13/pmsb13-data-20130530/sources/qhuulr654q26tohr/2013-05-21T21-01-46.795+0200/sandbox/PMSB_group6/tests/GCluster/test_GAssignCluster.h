#ifndef GINGER_TEST_GASSIGNCLUSTER_H_
#define GINGER_TEST_GASSIGNCLUSTER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GCluster/GAssignCluster.h>
#include <seqan/compareFiles.h>

using namespace seqan;
using namespace std;


SEQAN_DEFINE_TEST(test_GAssignCluster_dummy)
{
	GFastaRecord<DnaString> newRecord;
	newRecord.id="test";
	newRecord.seq="TTTCCAAA";
	GFastaRecord<DnaString> clusterRecord;
	newRecord.id="myClusterMaster";
	newRecord.seq="TTGGGGGGGGGGA";
	GFastaRecord<DnaString> newRecordRef;
	newRecordRef.id="test$myClusterMaster$";
	newRecordRef.seq="TTTCCAAA";
	
	SequenceStream outStream("../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/outFileGAssignCluster.fasta", SequenceStream::WRITE);
	SequenceStream outStreamRef("../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/outFileGAssignClusterRef.fasta", SequenceStream::WRITE);
	
	SEQAN_ASSERT(isGood(outStream));
	SEQAN_ASSERT(isGood(outStreamRef));
	SEQAN_ASSERT_NOT(writeRecord(outStreamRef, newRecordRef.id, newRecordRef.seq));
	SEQAN_ASSERT_NOT(GAssignCluster(outStream, newRecord, clusterRecord));
	
	FILE * pFile1;
	FILE * pFile2;
	pFile1 = fopen ("../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/outFileGAssignCluster.fasta" , "r");
	pFile2 = fopen ("../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/outFileGAssignClusterRef.fasta" , "r");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
}


#endif  // GINGER_TEST_GASSIGNCLUSTER_H_
