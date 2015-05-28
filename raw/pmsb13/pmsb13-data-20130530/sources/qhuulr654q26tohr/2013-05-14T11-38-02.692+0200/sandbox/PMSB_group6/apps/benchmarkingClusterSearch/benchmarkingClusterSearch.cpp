#include <seqan/GCluster.h>
#include <seqan/GSearch.h>
#include <seqan/GStructs/PerformanceSample.h>
#include <seqan/GStructs/MemorySample.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	String<GMatch<int> > matches;
	
	MemorySample ms("gesamt");
	
	const char * inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/rfam.10K.fasta";
	const char * inQueriesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_rfam.fasta";
	
	const char * outClustDataBasePath1 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.LA.fasta";
	const char * outMasterPath1 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.LA.fasta";
	const char * outMatchesPath1 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.LA.fasta";
	const char * scoreMode1 = "LOCAL_ALIGNMENT_MODE";
	
	const char * outClustDataBasePath2 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.GA.fasta";
	const char * outMasterPath2 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.GA.fasta";
	const char * outMatchesPath2 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.GA.fasta";
	const char * scoreMode2 = "GLOBAL_ALIGNMENT_MODE";
	
	const char * outClustDataBasePath3 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ.fasta";
	const char * outMasterPath3 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ.fasta";
	const char * outMatchesPath3 = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ.fasta";
	const char * scoreMode3 = "COMMON_QGRAM_MODE";
	
	const char * argvGCluster1[] = {"./bin/GCluster",inDataBasePath,outClustDataBasePath1,outMasterPath1,"-t","100","-m",scoreMode1};
	const char * argvGSearch1[] = {"./bin/GSearch",outClustDataBasePath3,inQueriesPath,outMatchesPath1,"-M",outMasterPath3,"-m",scoreMode1};
	const char * argvGCluster2[] = {"./bin/GCluster",inDataBasePath,outClustDataBasePath2,outMasterPath2,"-t","100","-m",scoreMode2};
	const char * argvGSearch2[] = {"./bin/GSearch",outClustDataBasePath3,inQueriesPath,outMatchesPath2,"-M",outMasterPath3,"-m",scoreMode2};
	const char * argvGCluster3[] = {"./bin/GCluster",inDataBasePath,outClustDataBasePath3,outMasterPath3,"-t","850","-m",scoreMode3};
	const char * argvGSearch3[] = {"./bin/GSearch",outClustDataBasePath3,inQueriesPath,outMatchesPath3,"-M",outMasterPath3,"-m",scoreMode3};
	int argcGSearch = 8;
	int argcGCluster = 8;
	
	
//	GCluster(performance, argcGCluster, argvGCluster1);
//	GCluster(performance, argcGCluster, argvGCluster2);
	GCluster(performance, argcGCluster, argvGCluster3);
	GSearch(matches, performance, argcGSearch, argvGSearch1);
	GSearch(matches, performance, argcGSearch, argvGSearch2);
	GSearch(matches, performance, argcGSearch, argvGSearch3);
	
	
	ms.computeMemory();
	for (int i=0; i<length(performance);++i) {
		performance[i].printToStdout();
	}
	ms.printToStdout();
	return 0;
}