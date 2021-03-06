#ifndef GINGER_GSEARCH_GSEARCH_BASE_H_
#define GINGER_GSEARCH_GSEARCH_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
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
	
	int res = readCommandlineParameters(clusteredDatabase, clusterRepresentanten, queries, outFile, scoreMode, threshold, argc, argv);
	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;
	
	bool MASTER_FILE_GIVEN = (clusterRepresentanten != toCString(""));
	
    SequenceStream myCluster(clusteredDatabase);
    SequenceStream myMaster(clusterRepresentanten);
    SequenceStream myQueries(queries);
    std::fstream myOutput(outFile, std::ios::binary | std::ios::out);
    if (checkStreams(myCluster, myMaster, myQueries, myOutput, MASTER_FILE_GIVEN))
        return 1;

    TScoringScheme scoringScheme(1, -2, -1);
	GFastaRecord<TSequence> queryRecord;
	GFastaRecord<TSequence> referenceRecord;
	GFastaRecord<TSequence> clusterRecord;

    TScore score;
	String<PerformanceSample> perfQueryAvgContainer;

	while (!readRecord(queryRecord.id, queryRecord.seq, myQueries)) {
		PerformanceSample perfQueryAvg("");
		open(myMaster, toCString(clusterRepresentanten));
		GScoreStorage<TScore> scoreStorage;
		
        while(MASTER_FILE_GIVEN && !readRecord(referenceRecord.id, referenceRecord.seq, myMaster)) {
			GScore(score, queryRecord, referenceRecord, scoringScheme, scoreMode);
			GEvaluateScore(scoreStorage, score, referenceRecord);
		}
        open(myCluster, toCString(clusteredDatabase));
		String<GMatch<TScore,TSequence> > matches;
		dataBaseLength=0;
        while(!readRecord(clusterRecord.id, clusterRecord.seq, myCluster)) {
			++dataBaseLength;
			if (GCheckClusterMaster(scoreStorage, clusterRecord, MASTER_FILE_GIVEN)) {
				GComputeMatch(matches, queryRecord, clusterRecord, scoringScheme, threshold, scoreMode);
            }
		}
		append(matchesCollection, matches);
        GPostProcessMatches(matches);
		GWriteMatches(myOutput, matches, queryRecord);
		
        perfQueryAvg.end();
		appendValue(perfQueryAvgContainer, perfQueryAvg);
    }
    PerformanceSample perfQueryAvg("Handle one query in GSearch average");
    perfQueryAvg.takeAverageValues(perfQueryAvgContainer);
	appendValue(performance, perfQueryAvg);
    perfSearchFunc.end();
	appendValue(performance, perfSearchFunc);
    return 0;
}

#endif  // GINGER_GSEARCH_GSEARCH_BASE_H_