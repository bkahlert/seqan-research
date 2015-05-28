//autor:Jakob

#include <seqan/GCluster.h>
#include <seqan/GSearch.h>
#include <seqan/GStructs/PerformanceSample.h>
#include <seqan/GStructs/MemorySample.h>
#include <seqan/computeSensSpec.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	
	const char * inDataBasePath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/pfam.100K.fasta";
	const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_pfam.included.fasta";
	
	const char * inDataBasePath2 = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/rfam.100K.fasta";
	const char * inQueriesPath2 = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/unmutable/queries_rfam.included.fasta";
	
	//const char * inDataBasePath = "../../seqan-trunk/sandbox/PMSB_group6/data/rfam_thresholdTest.fasta";
	//const char * inQueriesPath = "../../seqan-trunk/sandbox/PMSB_group6/data/rfam_thresholdTest_query.fasta";
	
	String<const char *> clustDatabases;
	appendValue(clustDatabases,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.pfam.CQ20.fasta");
	appendValue(clustDatabases,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.pfam.CQ40.fasta");
	appendValue(clustDatabases,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.pfam.CQ150.fasta"); 
	appendValue(clustDatabases,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.pfam.CQ500.fasta"); 
	
	String<const char *> clustDatabases2;
	appendValue(clustDatabases2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ20.fasta");
	appendValue(clustDatabases2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ40.fasta");
	appendValue(clustDatabases2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/clustered.rfam.CQ150.fasta"); 
	
	String<const char *> masterSequenceFiles;
	appendValue(masterSequenceFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.pfam.CQ20.fasta");
	appendValue(masterSequenceFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.pfam.CQ40.fasta");
	appendValue(masterSequenceFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.pfam.CQ150.fasta"); 
	appendValue(masterSequenceFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.pfam.CQ500.fasta"); 
	
	String<const char *> masterSequenceFiles2;
	appendValue(masterSequenceFiles2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ20.fasta");
	appendValue(masterSequenceFiles2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ40.fasta");
	appendValue(masterSequenceFiles2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/master.rfam.CQ150.fasta"); 
	
	String<int> clusterThresholds;
	appendValue(clusterThresholds, 20);
	appendValue(clusterThresholds, 40);
	appendValue(clusterThresholds, 150);
	appendValue(clusterThresholds, 500);
	
	String<const char *> matchFiles;
	appendValue(matchFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.pfam.CQ20.fasta");
	appendValue(matchFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.pfam.CQ40.fasta");
	appendValue(matchFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.pfam.CQ150.fasta");
	appendValue(matchFiles,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.pfam.CQ500.fasta");
	
	const char * matchFile_uncl = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.pfam.CQ_uncl.fasta";
	
	String<const char *> matchFiles2;
	appendValue(matchFiles2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ20.fasta");
	appendValue(matchFiles2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ40.fasta");
	appendValue(matchFiles2,"../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ150.fasta");
	
	const char * matchFile_uncl2 = "../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/matches.rfam.CQ_uncl.fasta";
	
	String<const char *> testRunNames;
	appendValue(testRunNames,"clustered GSearch t=20");
	appendValue(testRunNames,"clustered GSearch t=40");
	appendValue(testRunNames,"clustered GSearch t=150");
	appendValue(testRunNames,"clustered GSearch t=500");
	
	string scoreMode = "COMMON_QGRAM_MODE";
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	
	String<PerformanceSample> performance_uncl;
	String<GMatch<int, Peptide> > matchesSoll;
	int dataBaseLength;
	
	GSearch(matchesSoll, dataBaseLength, performance_uncl, inDataBasePath, "", inQueriesPath, matchFile_uncl, 500, scoreMode);
	performance_uncl[length(performance_uncl)-1].name="unclustered GSearch";
	double unclSpeed=performance_uncl[length(performance_uncl)-1].getTime();
	
	
	for (int i=0; i<length(performance_uncl);i+=1) {
		performance_uncl[i].printToStdout();
	}
	
	for (int i=0; i<length(clustDatabases)-1;++i){
		String<PerformanceSample> performance;
		String<GMatch<int, Peptide> > matchesIs;
		
		Peptide sequenceFormatExample;
		GCluster(performance,inDataBasePath,clustDatabases[i],masterSequenceFiles[i],sequenceFormatExample,clusterThresholds[i],0.5,scoreMode);
		
		GSearch(matchesIs, dataBaseLength, performance, clustDatabases[i], masterSequenceFiles[i], inQueriesPath, matchFiles[i], 500, scoreMode);
		performance[length(performance)-1].name=testRunNames[i];
		double clSpeed=performance[length(performance)-1].getTime();
		
		for (int i=0; i<length(performance);i+=1) {
			performance[i].printToStdout();
		}
		
		double sens;
		double spec;
		computeSensSpec(sens, spec, matchesSoll, matchesIs, dataBaseLength);
		cout << "Sensitivität:\t" << sens << endl;
		//cout << "Spezifität:\t" << spec << endl;
		cout << "Speedup:\t" << unclSpeed/clSpeed << endl << endl << endl;
	}
	
	cout << "===================================\n===============rfam================" << endl;
	//////////////////////////////////////////////////////////////////////////////////////////////
	
	String<PerformanceSample> performance_uncl2;
	String<GMatch<int, Rna5String> > matchesSoll2;
	int dataBaseLength2;
	
	GSearch(matchesSoll2, dataBaseLength2, performance_uncl2, inDataBasePath2, "", inQueriesPath2, matchFile_uncl2, 400, scoreMode);
	performance_uncl2[length(performance_uncl2)-1].name="unclustered GSearch";
	unclSpeed=performance_uncl2[length(performance_uncl2)-1].getTime();
	
	
	for (int i=0; i<length(performance_uncl2);i+=1) {
		performance_uncl2[i].printToStdout();
	}
	
	for (int i=0; i<length(clustDatabases2)-1;++i){
		String<PerformanceSample> performance2;
		String<GMatch<int, Rna5String> > matchesIs2;
		
		Rna5String sequenceFormatExample2;
		GCluster(performance2,inDataBasePath2,clustDatabases2[i],masterSequenceFiles2[i],sequenceFormatExample2,clusterThresholds[i],0.5,scoreMode);
		
		GSearch(matchesIs2, dataBaseLength2, performance2, clustDatabases2[i], masterSequenceFiles2[i], inQueriesPath2, matchFiles2[i], 400, scoreMode);
		performance2[length(performance2)-1].name=testRunNames[i];
		double clSpeed=performance2[length(performance2)-1].getTime();
		
		for (int i=0; i<length(performance2);i+=1) {
			performance2[i].printToStdout();
		}
		
		double sens2;
		double spec2;
		computeSensSpec(sens2, spec2, matchesSoll2, matchesIs2, dataBaseLength2);
		cout << "Sensitivität:\t" << sens2 << endl;
		//cout << "Spezifität:\t" << spec << endl;
		cout << "Speedup:\t" << unclSpeed/clSpeed << endl << endl << endl;
	}
	
	return 0;
}