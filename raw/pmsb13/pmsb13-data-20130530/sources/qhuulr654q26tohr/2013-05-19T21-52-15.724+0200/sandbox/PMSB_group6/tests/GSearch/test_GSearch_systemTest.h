#ifndef GINGER_TEST_GSEARCH_SYSTEMTEST_H_
#define GINGER_TEST_GSEARCH_SYSTEMTEST_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>
#include <seqan/compareFiles.h>

using namespace seqan;
using namespace std;

// empty files & LOCAL_ALIGNMENT_MODE.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_emptyFiles_LAMode)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_empty_reference.fasta";
	const char * inMasterReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_empty_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_empty.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_empty_LocalAlignmentMode.txt";
	const char * outMatchesReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/matches_empty_reference.txt";
	string scoreMode = "LOCAL_ALIGNMENT_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Dna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, inMasterReferencePath, inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
	
	FILE * pFile1;
	FILE * pFile2;
	pFile1 = fopen (outMatchesPath , "r");
	if (pFile1 == NULL)
		SEQAN_FAIL("Error: Could not open output path");
	pFile2 = fopen (outMatchesReferencePath , "r");
	if (pFile2 == NULL)
		SEQAN_FAIL("Error: Could not open output reference path");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
}

// empty files & COMMON_QGRAM_MODE.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_emptyFiles_CQMode)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_empty_reference.fasta";
	const char * inMasterReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_empty_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_empty.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_empty_CommonQgramMode.txt";
	const char * outMatchesReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/matches_empty_reference.txt";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Dna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, inMasterReferencePath, inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
	
	FILE * pFile1;
	FILE * pFile2;
	pFile1 = fopen (outMatchesPath , "r");
	if (pFile1 == NULL)
		SEQAN_FAIL("Error: Could not open output path");
	pFile2 = fopen (outMatchesReferencePath , "r");
	if (pFile2 == NULL)
		SEQAN_FAIL("Error: Could not open output reference path");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
}

// empty files & COMMON_QGRAM_MODE & unclustered.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_emptyFiles_CQMode_unclustered)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_empty_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_empty.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_empty_CommonQgramMode_unclustered.txt";
	const char * outMatchesReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/matches_empty_reference.txt";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Dna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, "", inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
	
	FILE * pFile1;
	FILE * pFile2;
	pFile1 = fopen (outMatchesPath , "r");
	if (pFile1 == NULL)
		SEQAN_FAIL("Error: Could not open output path");
	pFile2 = fopen (outMatchesReferencePath , "r");
	if (pFile2 == NULL)
		SEQAN_FAIL("Error: Could not open output reference path");
	SEQAN_ASSERT_EQ(compareFiles(pFile1,pFile2),0);
}


// toyExample files & LOCAL_ALIGNMENT_MODE.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_toyFiles_LAMode)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_toy_reference.fasta";
	const char * inMasterReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_toy_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_toy.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_toy_LocalAlignmentMode.txt";
	string scoreMode = "LOCAL_ALIGNMENT_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Dna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, inMasterReferencePath, inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
}

// toyExample files & COMMON_QGRAM_MODE.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_toyFiles_CQMode)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_toy_reference.fasta";
	const char * inMasterReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_toy_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_toy.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_toy_CommonQgramMode.txt";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Dna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, inMasterReferencePath, inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
}

// toyExample files & COMMON_QGRAM_MODE & unclustered.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_toyFiles_CQMode_unclustered)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_toy_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_toy.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_toy_CommonQgramMode_unclustered.txt";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Dna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, inMasterReferencePath, inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
}


// rfam_short files & LOCAL_ALIGNMENT_MODE.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_rfamFiles_LAMode)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_rfam_reference.fasta";
	const char * inMasterReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_rfam_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_rfam.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_rfam_LocalAlignmentMode.txt";
	string scoreMode = "LOCAL_ALIGNMENT_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Rna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, inMasterReferencePath, inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
}

// rfam_short files & LOCAL_ALIGNMENT_MODE & unclustered.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_rfamFiles_LAMode_unclustered)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_rfam_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_rfam.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_rfam_LocalAlignmentMode_unclustered.txt";
	string scoreMode = "LOCAL_ALIGNMENT_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Rna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, inMasterReferencePath, inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
}

// rfam_short files & COMMON_QGRAM_MODE.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_rfamFiles_CQMode)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/clusteredDatabase_rfam_reference.fasta";
	const char * inMasterReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/master_rfam_reference.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_rfam.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches_rfam_CommonQgramMode.txt";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Rna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, inMasterReferencePath, inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
}

// unbekannter Speicherzugriffsfehler.
SEQAN_DEFINE_TEST(test_GSearch_systemTest_speicherzugriffsfehler_0517)
{
	const char * inClustDataBaseReferencePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/speicherzugriffsfehler_0517.dataBase.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/speicherzugriffsfehler_0517.query.fasta";
	const char * outMatchesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/speicherzugriffsfehler_0517.matches.txt";
	string scoreMode = "COMMON_QGRAM_MODE";
	
	String<PerformanceSample> performance;
	String<GMatch<int, Rna5String> > matches;
	int dataBaseLength;
	
	int res = GSearch(matches, dataBaseLength, performance, inClustDataBaseReferencePath, "", inQueriesPath, outMatchesPath,1,scoreMode);
	
	SEQAN_ASSERT_EQ(res,0);
}

#endif  // GINGER_TEST_GSEARCH_SYSTEMTEST_H_
