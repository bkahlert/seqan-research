#include <seqan/clusterSearch.h>
#include <seqan/structs/PerformanceSample.h>

using namespace seqan;
using namespace std;

int main(int argcm, char const ** argmv) {
	
	const char * outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/clusteredDatabase.fasta";
	const char * outClustDataBaseReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/clusteredDatabase_reference.fasta";
	const char * outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/master.fasta";
	const char * outMasterReferencePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/master_reference.fasta";
	
	const char * argv[] = {"./bin/clusterSearch","/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/emptyTest/empty.fasta",outClustDataBasePath,outMasterPath,"-t","4"};
	int argc = 6;
	String<PerformanceSample> performance;
	cout << argc << endl;
	for (int i=0;i<argc;++i)
		cout << "-> " << argv[i] << endl;
	cout << "----------" << endl;
	int res = clusterSearch(performance, argc, argv);
	
	/*
	cout << argc << endl;
	for (int i=0;i<argc;++i)
		cout << argv[i] << endl;
	String<PerformanceSample> performance;
	
	int res = clusterSearch(performance, argc, argv);
	return res;*/
}