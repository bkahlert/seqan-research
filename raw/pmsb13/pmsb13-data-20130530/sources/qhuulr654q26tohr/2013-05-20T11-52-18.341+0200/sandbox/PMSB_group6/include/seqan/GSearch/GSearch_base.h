//Autor:Jakob

#ifndef GINGER_GSEARCH_GSEARCH_BASE_H_
#define GINGER_GSEARCH_GSEARCH_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fstream>
#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/GMatch.h>
#include <seqan/GStructs/GScoreStorage.h>
#include <seqan/GStructs/PerformanceSample.h>
#include <seqan/GStructs/MemorySample.h>
#include <seqan/GSearch/GScore.h>
#include <seqan/GSearch/GEvaluateScore.h>
#include <seqan/GSearch/GParseMasterId.h>
#include <seqan/GSearch/GCheckClusterMaster.h>
#include <seqan/GSearch/GComputeMatch.h>
#include <seqan/GSearch/GPostProcessMatches.h>
#include <seqan/GSearch/GWriteMatches.h>

using namespace seqan;
using namespace std;


int checkStreams(SequenceStream &myCluster, SequenceStream &myMaster, SequenceStream &myQueries, fstream &myOutput, bool MASTER_FILE_GIVEN) {
    if (!isGood(myCluster)) {
        std::cerr << "ERROR: Could not open the clustered database file in function GSearch.\n";
        return 1;
    }

    if (!isGood(myMaster) && MASTER_FILE_GIVEN) {
		std::cerr << "ERROR: Could not open the master sequences file in function GSearch.\n";
        return 1;
    }
    if (!isGood(myQueries)) {
		std::cerr << "ERROR: Could not open the queries file in function GSearch.\n";
        return 1;
    }
    if (!myOutput.good()) {
		std::cerr << "ERROR: Could not open the output file in function GSearch.\n";
        return 1;
    }
    return 0;
}


template <typename TScore, typename TSequence>
int GSearch(String<GMatch<TScore,TSequence> > & matchesCollection, int & dataBaseLength, String<PerformanceSample> & performance,
			const char * clusteredDatabase, const char * clusterRepresentanten, const char * queries, const char * outFile, TScore threshold, string scoreMode) {
	
	PerformanceSample perfSearchFunc("GSearch Funktion");
    typedef Score<int> TScoringScheme;
	
	bool MASTER_FILE_GIVEN = (clusterRepresentanten != toCString(""));
	
    SequenceStream myCluster(clusteredDatabase);
    SequenceStream myMaster(clusterRepresentanten);
    SequenceStream myQueries(queries);
    std::fstream myOutput(outFile, std::ios::binary | std::ios::out);
    if (checkStreams(myCluster, myMaster, myQueries, myOutput, MASTER_FILE_GIVEN))
        return 1;

	//cout << "open streams completed" << endl;
    TScoringScheme scoringScheme(1, -2, -1);
	GFastaRecord<TSequence> queryRecord;
	GFastaRecord<TSequence> referenceRecord;
	GFastaRecord<TSequence> clusterRecord;

    TScore score;
	String<PerformanceSample> perfQueryAvgContainer;

	while (!readRecord(queryRecord.id, queryRecord.seq, myQueries)) {
		//cout << "read record completed" << endl;
		PerformanceSample perfQueryAvg("");
		open(myMaster, toCString(clusterRepresentanten));
		GScoreStorage<TScore> scoreStorage;
		
        while(MASTER_FILE_GIVEN && !readRecord(referenceRecord.id, referenceRecord.seq, myMaster)) {
			GScore(score, queryRecord, referenceRecord, scoringScheme, scoreMode);
			GEvaluateScore(scoreStorage, score, referenceRecord);
		}
		//cout << "find master completed" << endl;
		//cout << scoreStorage.maxId << endl;
		//cout << scoreStorage.maxScore << endl;
        open(myCluster, toCString(clusteredDatabase));
		String<GMatch<TScore,TSequence> > matches;
		dataBaseLength=0;
		while(!readRecord(clusterRecord.id, clusterRecord.seq, myCluster)) {
			//cout << "read record: " << clusterRecord.seq << endl;
			++dataBaseLength;
			if (GCheckClusterMaster(scoreStorage, clusterRecord, MASTER_FILE_GIVEN)) {
				//cout << "compute match ..." << endl;
				GComputeMatch(matches, queryRecord, clusterRecord, scoringScheme, threshold, scoreMode);
            }
		}
		//cout << "find matches completed" << endl;
		append(matchesCollection, matches);
        GPostProcessMatches(matches);
		GWriteMatches(myOutput, matches, queryRecord);
		
		//cout << "Searched one query" << endl;
        perfQueryAvg.end();
		appendValue(perfQueryAvgContainer, perfQueryAvg);
	}
	cout << "\tSearched "<<length(perfQueryAvgContainer)<<" queries and found "<<length(matches)<<" matches."
    PerformanceSample perfQueryAvg("Handle one query in GSearch average");
    perfQueryAvg.takeAverageValues(perfQueryAvgContainer);
	appendValue(performance, perfQueryAvg);
    perfSearchFunc.end();
	appendValue(performance, perfSearchFunc);
    return 0;
}

#endif  // GINGER_GSEARCH_GSEARCH_BASE_H_