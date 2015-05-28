#ifndef SANDBOX_PMSB_GROUP6_TEST_CLUSTERSEARCH_H_
#define SANDBOX_PMSB_GROUP6_TEST_CLUSTERSEARCH_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/clusterSearch.h>
#include <seqan/structs/PerformanceSample.h>
#include <stdio.h>

using namespace seqan;
using namespace std;

int compareFiles(FILE* f1, FILE* f2) {
	int N = 10000;
	char buf1[N];
	char buf2[N];
	
	do {
		size_t r1 = fread(buf1, 1, N, f1);
		size_t r2 = fread(buf2, 1, N, f2);
		
		if (r1 != r2 ||
			memcmp(buf1, buf2, r1)) {
			return 1;  // Files are not equal
			}
	} while (!feof(f1) && !feof(f2));
	
	return !(feof(f1) && feof(f2));
}

// A test for strings.
SEQAN_DEFINE_TEST(test_clusterSearch_emptyFiles)
{
	const char * outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/clusteredDatabase.fasta";
	const char * outClustDataBaseReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/clusteredDatabase_reference.fasta";
	const char * outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/master.fasta";
	const char * outMasterReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/master_reference.fasta";
	
	const char * argv[] = {"./bin/clusterSearch","/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/empty.fasta",outClustDataBasePath,outMasterPath,"-t","4"};
	int argc = 6;
	String<PerformanceSample> performance;
	int res = clusterSearch(performance, argc, argv);
	
	
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

SEQAN_DEFINE_TEST(test_clusterSearch_benchmarking_emptyFiles)
{	
	SEQAN_FAIL("not yet implemented");
}

#endif  // SANDBOX_PMSB_GROUP6_TESTS_CLUSTERSEARCH_TEST_CLUSTERSEARCH_H_
