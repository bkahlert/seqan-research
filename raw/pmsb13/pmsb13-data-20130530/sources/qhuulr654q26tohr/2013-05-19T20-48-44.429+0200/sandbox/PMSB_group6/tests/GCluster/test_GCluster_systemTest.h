#ifndef SANDBOX_PMSB_GROUP6_TEST_CLUSTER_SYSTEMTEST_H_
#define SANDBOX_PMSB_GROUP6_TEST_CLUSTER_SYSTEMTEST_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GCluster.h>
#include <seqan/compareFiles.h>

using namespace seqan;
using namespace std;

// empty files & LOCAL_ALIGNMENT_MODE.
SEQAN_DEFINE_TEST(test_GCluster_systemTest_emptyFiles_LAMode)
{
	string inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/empty.fasta";
	string outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clusteredDatabase_empty_LocalAlignmentMode.fasta";
	string outClustDataBaseReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_empty_reference.fasta";
	string outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master_empty_LocalAlignmentMode.fasta";
	string outMasterReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_empty_reference.fasta";
	string scoreMode = "LOCAL_ALIGNMENT_MODE";
	
	String<PerformanceSample> performance;
	DnaString sequenceFormatExample;
	int res = GCluster(performance, inDataBasePath, outClustDataBasePath, outMasterPath, sequenceFormatExample, 4, 0.5, scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
	
	FILE * pFile1;
	FILE * pFile2;
	pFile1 = fopen (outClustDataBasePath.c_str() , "r");
	pFile2 = fopen (outClustDataBaseReferencePath.c_str() , "r");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
	
	pFile1 = fopen (outMasterPath.c_str() , "r");
	pFile2 = fopen (outMasterReferencePath.c_str() , "r");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
}

// empty files & COMMON_QGRAM_MODE.
SEQAN_DEFINE_TEST(test_GCluster_systemTest_emptyFiles_CQMode)
{
	string inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/empty.fasta";
	string outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clusteredDatabase_empty.fasta";
	string outClustDataBaseReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_empty_reference.fasta";
	string outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master_empty.fasta";
	string outMasterReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_empty_reference.fasta";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	int res = GCluster(performance, inDataBasePath, outClustDataBasePath, outMasterPath, sequenceFormatExample, 800, 0.5, scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
	
	FILE * pFile1;
	FILE * pFile2;
	pFile1 = fopen (outClustDataBasePath , "r");
	pFile2 = fopen (outClustDataBaseReferencePath , "r");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
	
	pFile1 = fopen (outMasterPath , "r");
	pFile2 = fopen (outMasterReferencePath , "r");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
}



// toy files & LOCAL_ALIGNMENT_MODE.
SEQAN_DEFINE_TEST(test_GCluster_systemTest_toyFiles_LAMode)
{
	
	string inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/toyDatabase.fasta";
	string outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clusteredDatabase_toy.fasta";
	string outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master_toy.fasta";
	string scoreMode = "LOCAL_ALIGNMENT_MODE";
	
	String<PerformanceSample> performance;
	DnaString sequenceFormatExample;
	int res = GCluster(performance, inDataBasePath, outClustDataBasePath, outMasterPath, sequenceFormatExample, 4, 0.5, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
}

// toy files & COMMON_QGRAM_MODE.
SEQAN_DEFINE_TEST(test_GCluster_systemTest_toyFiles_CQMode)
{
	
	string inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/toyDatabase.fasta";
	string outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clusteredDatabase_toy.fasta";
	string outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master_toy.fasta";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	String<PerformanceSample> performance;
	DnaString sequenceFormatExample;
	int res = GCluster(performance, inDataBasePath, outClustDataBasePath, outMasterPath, sequenceFormatExample, 850, 0.5, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
}

// rfam_short files & LOCAL_ALIGNMENT_MODE.
SEQAN_DEFINE_TEST(test_GCluster_systemTest_rfamFiles_LAMode)
{
	string inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/rfam.10K.fasta";
	string outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clusteredDatabase_rfam.fasta";
	string outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master_rfam.fasta";
	string scoreMode = "LOCAL_ALIGNMENT_MODE";
	
	String<PerformanceSample> performance;
	DnaString sequenceFormatExample;
	int res = GCluster(performance, inDataBasePath, outClustDataBasePath, outMasterPath, sequenceFormatExample, 100, 0.5, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
}

// rfam_short files & COMMON_QGRAM_MODE.
SEQAN_DEFINE_TEST(test_GCluster_systemTest_rfamFiles_CQMode)
{
	string inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/rfam.10K.fasta";
	string outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clusteredDatabase_rfam.fasta";
	string outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master_rfam.fasta";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	String<PerformanceSample> performance;
	DnaString sequenceFormatExample;
	int res = GCluster(performance, inDataBasePath, outClustDataBasePath, outMasterPath, sequenceFormatExample, 850, 0.5, scoreMode);
	SEQAN_ASSERT_EQ(res,0);
}


#endif  // SANDBOX_PMSB_GROUP6_TEST_CLUSTER_SYSTEMTEST_H_
