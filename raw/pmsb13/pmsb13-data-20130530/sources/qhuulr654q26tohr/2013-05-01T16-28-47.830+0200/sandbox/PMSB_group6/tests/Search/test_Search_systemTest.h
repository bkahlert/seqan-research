#ifndef SANDBOX_PMSB_GROUP6_TEST_SEARCH_SYSTEMTEST_H_
#define SANDBOX_PMSB_GROUP6_TEST_SEARCH_SYSTEMTEST_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/Search.h>

using namespace seqan;
using namespace std;

// returns true if two files are exactly the same
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

// empty files.
SEQAN_DEFINE_TEST(test_Search_systemTest_emptyFiles)
{
	const char * inClustDataBaseReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_empty_reference.fasta";
	const char * inMasterReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_empty_reference.fasta";
	const char * inQueriesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_empty.fasta";
	const char * outMatchesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_empty.txt";
	const char * outMatchesReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/matches_empty_reference.txt";
	
	const char * argv[] = {"./bin/Search",inClustDataBaseReferencePath,inMasterReferencePath,inQueriesPath,outMatchesPath};
	int argc = 5;
	String<PerformanceSample> performance;
	int res = Search(performance, argc, argv);
	
	
	SEQAN_ASSERT_EQ(res,0);
	
	FILE * pFile1;
	FILE * pFile2;
	pFile1 = fopen (outMatchesPath , "r");
	if (pFile1 == NULL){
		cerr << "Error: Could not open " << outMatchesPath << endl;
		SEQAN_FAIL;}
	pFile2 = fopen (outMatchesReferencePath , "r");
	if (pFile2 == NULL){
		cerr << "Error: Could not open " << outMatchesPath << endl;
		SEQAN_FAIL;}
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
}


// toyExample files.
SEQAN_DEFINE_TEST(test_Search_systemTest_toyFiles)
{
	const char * inClustDataBaseReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_toy_reference.fasta";
	const char * inMasterReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_toy_reference.fasta";
	const char * inQueriesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_toy.fasta";
	const char * outMatchesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_toy.txt";
	const char * outMatchesReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_toy_reference.txt";
	
	const char * argv[] = {"./bin/Search",inClustDataBaseReferencePath,inMasterReferencePath,inQueriesPath,outMatchesPath};
	int argc = 5;
	String<PerformanceSample> performance;
	int res = Search(performance, argc, argv);
	
	
	SEQAN_ASSERT_EQ(res,0);
	
	FILE * pFile1;
	FILE * pFile2;
	pFile1 = fopen (outMatchesPath , "r");
	pFile2 = fopen (outMatchesReferencePath , "r");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
}

#endif  // SANDBOX_PMSB_GROUP6_TEST_SEARCH_SYSTEMTEST_H_
