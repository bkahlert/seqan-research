#ifndef SANDBOX_PMSB_GROUP6_TEST_CLUSTER_SYSTEMTEST_H_
#define SANDBOX_PMSB_GROUP6_TEST_CLUSTER_SYSTEMTEST_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/cluster.h>
#include <seqan/compareFiles.h>

using namespace seqan;
using namespace std;

// empty files.
SEQAN_DEFINE_TEST(test_cluster_systemTest_emptyFiles)
{
	const char * outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clusteredDatabase_empty.fasta";
	const char * outClustDataBaseReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_empty_reference.fasta";
	const char * outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master_empty.fasta";
	const char * outMasterReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_empty_reference.fasta";
	
	const char * argv[] = {"./bin/cluster","/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/empty.fasta",outClustDataBasePath,outMasterPath,"-t","4"};
	int argc = 6;
	String<PerformanceSample> performance;
	int res = cluster(performance, argc, argv);
	
	
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


// toyExample files.
SEQAN_DEFINE_TEST(test_cluster_systemTest_toyFiles)
{
	const char * inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/toyDatabase.fasta";
	const char * outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clusteredDatabase_toy.fasta";
	const char * outClustDataBaseReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_toy_reference.fasta";
	const char * outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master_toy.fasta";
	const char * outMasterReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_toy_reference.fasta";
	
	const char * argv[] = {"./bin/cluster",inDataBasePath,outClustDataBasePath,outMasterPath,"-t","4"};
	int argc = 6;
	String<PerformanceSample> performance;
	int res = cluster(performance, argc, argv);
	
	
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

#endif  // SANDBOX_PMSB_GROUP6_TEST_CLUSTER_SYSTEMTEST_H_
