#include <seqan/cluster.h>
#include <seqan/Search.h>
#include <seqan/structs/PerformanceSample.h>
#include <seqan/structs/MemorySample.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	
	MemorySample ms("gesamt");
	
	const char * inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/data/rfam.1M.fasta";
	const char * outClustDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/data/clustered.rfam.fasta";
	const char * outMasterPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/data/master.rfam.fasta";
	
	const char * argv1[] = {"./bin/cluster",inDataBasePath,outClustDataBasePath,outMasterPath,"-t","4000"};
	int argc1 = 6;
	
	const char * inClustDataBaseReferencePath = outClustDataBasePath;
	const char * inMasterReferencePath = outMasterPath;
	const char * inQueriesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/data/rfam.queries.fasta";
	const char * outMatchesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/data/matches.rfam.txt";
	
	const char * argv2[] = {"./bin/Search",inClustDataBaseReferencePath,inMasterReferencePath,inQueriesPath,outMatchesPath};
	int argc2 = 5;
	
	int res1 = cluster(performance, argc1, argv1);
	int res2 = Search(performance, argc2, argv2);
	ms.computeMemory();
	
	for (int i=0; i<length(performance);++i) {
		performance[i].printToStdout();
	}
	ms.printToStdout();
	return 0;
}