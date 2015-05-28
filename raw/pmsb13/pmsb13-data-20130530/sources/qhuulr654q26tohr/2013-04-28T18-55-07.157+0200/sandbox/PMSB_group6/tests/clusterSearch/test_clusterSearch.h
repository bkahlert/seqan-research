#ifndef SANDBOX_PMSB_GROUP6_TEST_CLUSTERSEARCH_H_
#define SANDBOX_PMSB_GROUP6_TEST_CLUSTERSEARCH_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/clusterSearch.h>
#include <seqan/PerformanceSample.h

using namespace seqan;
using namespace std;

// A test for strings.
SEQAN_DEFINE_TEST(test_clusterSearch_emptyFiles)
{
	char ** argv = {"clusterSearch","../autoTestData/empty.fasta","../autoTestData/clusteredDatabase.fasta","../autoTestData/master.fasta","-t","4"}
	int argc = length(argv);
	String<PerformanceSample> performance;
	int res = clusterSearch(performance, argc, argv);
}

SEQAN_DEFINE_TEST(test_clusterSearch_benchmarking_emptyFiles)
{	
	SEQAN_FAIL("not yet implemented");
}

#endif  // SANDBOX_PMSB_GROUP6_TESTS_CLUSTERSEARCH_TEST_CLUSTERSEARCH_H_
