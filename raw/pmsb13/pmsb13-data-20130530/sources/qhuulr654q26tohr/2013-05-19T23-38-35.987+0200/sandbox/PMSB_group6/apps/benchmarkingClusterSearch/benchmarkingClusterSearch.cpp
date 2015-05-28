#include <seqan/GCluster.h>
#include <seqan/GSearch.h>
#include <seqan/GStructs/PerformanceSample.h>
#include <seqan/GStructs/MemorySample.h>
#include <seqan/computeSensSpec.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	
	const char * inDataBasePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/rfam.100K.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_rfam.short.fasta";
	
	//const char * inDataBasePath = "../../seqan-trunk/sandbox/PMSB_group6/data/rfam_thresholdTest.fasta";
	//const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/data/rfam_thresholdTest_query.fasta";
	
	String<const char *> clustDatabases;
	appendValue(clustDatabases,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ1.fasta");
	appendValue(clustDatabases,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ2.fasta");
	appendValue(clustDatabases,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ3.fasta"); 
	 
	String<const char *> masterSequenceFiles;
	appendValue(masterSequenceFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ1.fasta");
	appendValue(masterSequenceFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ2.fasta");
	appendValue(masterSequenceFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ3.fasta"); 
	
	String<int> clusterThresholds;
	appendValue(clusterThresholds, 20);
	appendValue(clusterThresholds, 100);
	appendValue(clusterThresholds, 200);
	
	String<const char *> matchFiles;
	appendValue(matchFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ1.fasta");
	appendValue(matchFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ2.fasta");
	appendValue(matchFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ3.fasta");
	
	const char * matchFile_uncl = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ_uncl.fasta";
	
	String<const char *> testRunNames;
	appendValue(testRunNames,"clustered GSearch");
	appendValue(testRunNames,"matches.rfam.CQ2.fasta");
	appendValue(testRunNames,"matches.rfam.CQ3.fasta");
	
	string scoreMode = "COMMON_QGRAM_MODE";
	
	
	
	String<PerformanceSample> performance_uncl;
	String<GMatch<int, Rna5String> > matchesSoll;
	int dataBaseLength;
	
	GSearch(matchesSoll, dataBaseLength, performance_uncl, inDataBasePath, "", inQueriesPath, matchFile_uncl, 100, scoreMode);
	performance_uncl[length(performance_uncl)-1].name="unclustered GSearch";
	double unclSpeed=performance_uncl[length(performance_uncl)-1].getTime();
	
	
	for (int i=0; i<length(performance_uncl);++i) {
		performance_uncl[i].printToStdout();
	}
	
	for (int i=0; i<length(clustDatabases)-2;++i){
		String<PerformanceSample> performance;
		String<GMatch<int, Rna5String> > matchesIs;
		
		Rna5String sequenceFormatExample;
		GCluster(performance,inDataBasePath,clustDatabases[i],masterSequenceFiles[i],sequenceFormatExample,clusterThresholds[i],0.5,scoreMode);
		
		GSearch(matchesIs, dataBaseLength, performance, clustDatabases[i], masterSequenceFiles[i], inQueriesPath, matchFiles[i], 100, scoreMode);
		performance[length(performance)-1].name=testRunNames[i];
		double clSpeed=performance[length(performance)-1].getTime();
		
		for (int i=0; i<length(performance);++i) {
			performance[i].printToStdout();
		}
		
		double sens;
		double spec;
		computeSensSpec(sens, spec, matchesSoll, matchesIs, dataBaseLength);
		cout << "Sensitivität:\t" << sens << endl;
		//cout << "Spezifität:\t" << spec << endl;
		cout << "Speedup:\t" << unclSpeed/clSpeed << endl << endl << endl;
	}
	
	return 0;
}