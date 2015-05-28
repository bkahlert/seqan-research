#include <seqan/GCluster.h>
#include <seqan/GSearch.h>
#include <seqan/GStructs/PerformanceSample.h>
#include <seqan/GStructs/MemorySample.h>
#include <seqan/computeSensSpec.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	
	//const char * inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/rfam.100K.fasta";
	//const char * inQueriesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_rfam.short.fasta";
	
	const char * inDataBasePath = "/home/development/seqan-trunk/sandbox/PMSB_group6/data/rfam_thresholdTest.fasta";
	const char * inQueriesPath = "/home/development/seqan-trunk/sandbox/PMSB_group6/data/error_query_thresholdTest.fasta";
	
	String<const char *> clustDatabases;
	appendValue(clustDatabases,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ1.fasta");
	appendValue(clustDatabases,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ2.fasta");
	appendValue(clustDatabases,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ3.fasta"); 
	
	String<const char *> masterSequenceFiles;
	appendValue(masterSequenceFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ1.fasta");
	appendValue(masterSequenceFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ2.fasta");
	appendValue(masterSequenceFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ3.fasta"); 
	
	String<const char *> clusterThresholds;
	appendValue(clusterThresholds, "20");
	appendValue(clusterThresholds, "100");
	appendValue(clusterThresholds, "200");
	
	String<const char *> matchFiles;
	appendValue(matchFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ1.1.fasta");
	appendValue(matchFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ1.2.fasta");
	appendValue(matchFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ2.1.fasta");
	appendValue(matchFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ2.2.fasta");
	appendValue(matchFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ3.1.fasta");
	appendValue(matchFiles,"/home/development/seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ3.2.fasta");
	
	String<const char *> testRunNames;
	appendValue(testRunNames,"matches.rfam.CQ1.1.fasta");
	appendValue(testRunNames,"matches.rfam.CQ1.2.fasta");
	appendValue(testRunNames,"matches.rfam.CQ2.1.fasta");
	appendValue(testRunNames,"matches.rfam.CQ2.2.fasta");
	appendValue(testRunNames,"matches.rfam.CQ3.1.fasta");
	appendValue(testRunNames,"matches.rfam.CQ3.2.fasta");
	
	const char * scoreMode = "COMMON_QGRAM_MODE";
	
	for (int i=0; i<length(clustDatabases);++i){
		String<PerformanceSample> performance;
		String<GMatch<int> > matchesIs;
		String<GMatch<int> > matchesSoll;
		int dataBaseLength;
		
		const char * argvGCluster[] = {"./bin/GCluster",inDataBasePath,clustDatabases[i],masterSequenceFiles[i],"-t",clusterThresholds[i],"-m",scoreMode};
		int argcGCluster = 8;
		GCluster(performance, argcGCluster, argvGCluster);
		
		const char * argvGSearch[] = {"./bin/GSearch",clustDatabases[i],inQueriesPath,matchFiles[(i*2)],"-M",masterSequenceFiles[i],"-m",scoreMode};
		int argcGSearch = 8;
		GSearch(matchesIs, dataBaseLength, performance, argcGSearch, argvGSearch);
		performance[length(performance)-1].name=testRunNames[(i*2)];
		
		const char * argvGSearch2[] = {"./bin/GSearch",clustDatabases[i],inQueriesPath,matchFiles[(i*2)+1],"-m",scoreMode};
		GSearch(matchesSoll, dataBaseLength, performance, argcGSearch-2, argvGSearch2);
		performance[length(performance)-1].name=testRunNames[(i*2)+1];
		
		
		for (int i=0; i<length(performance);++i) {
			performance[i].printToStdout();
		}
		
		double sens;
		double spec;
		computeSensSpec(sens, spec, matchesSoll, matchesIs, dataBaseLength);
		cout << "Sensitivität:\t" << sens << endl;
		cout << "Spezifität:\t" << spec << endl << endl << endl;
	}
	
	
	return 0;
}