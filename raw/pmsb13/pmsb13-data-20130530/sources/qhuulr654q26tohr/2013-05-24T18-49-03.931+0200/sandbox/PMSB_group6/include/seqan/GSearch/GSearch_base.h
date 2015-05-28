//Autor:Hannes

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
#include <seqan/GSearch/GScorePeptideQGram.h>
#include <seqan/GSearch/GEvaluateScore.h>
#include <seqan/GSearch/GParseMasterId.h>
#include <seqan/GSearch/GCheckClusterMaster.h>
#include <seqan/GSearch/GComputeMatch.h>
#include <seqan/GSearch/GPostProcessMatches.h>
#include <seqan/GSearch/GWriteMatches.h>

using namespace seqan;
using namespace std;

/**
 * \brief Checks if the given streams are good.
 */
int checkStreams(SequenceStream &myCluster/**<[in]path of the cluster file*/, SequenceStream &myMaster/**<[in]path of the master file*/, SequenceStream &myQueries/**<[in]path of the query file*/, fstream &myOutput/**<[in]path of the output file*/, bool MASTER_FILE_GIVEN/**<[in]specifies if thesearch algorithm runs on a clustered or unclustered database*/) {
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
/**
 * \brief Searches for a sequence similar to a given query in a clustered database.
 * 
 * This is the base file of the GSearch algorithm. It includes all relevant functions and executes them.
 * Therefore the algorithm searches every query, and for every query it searches the ClusterMaster file to determine the Mastersequence with the highest similarity.
 * For every sequence in the given cluster it determines if the match score is higher then the threshold score and, if that is given, builds a global alignment with the sequence in the cluster and the query and adds the result to the matches.  
 * The results get sorted after that.
 */
int GSearch(String<GMatch<TScore,TSequence> > & matchesCollection/**<[out]container of the found matches*/, int & dataBaseLength/**<[out]container for length of the database*/, String<PerformanceSample> & performance/**<[out]Struct in which the performance data is written*/,
			const char * clusteredDatabase/**<[in]the database to be searched in*/, const char * clusterRepresentanten/**<[in]the file with all the clusterMaster*/, const char * queries/**<[in]the queries to be searched*/, const char * outFile/**<[in]the path of the file to be written in*/, TScore threshold/**<[in]used threshold*/, string scoreMode/**<[in]used Scoring Mode*/) {
	
	PerformanceSample perfSearchFunc("GSearch Funktion");
    typedef Score<int> TScoringScheme;
	
	CharString clusterRepresentantenString=clusterRepresentanten;
	bool MASTER_FILE_GIVEN = (clusterRepresentantenString != "NOFILE" && clusterRepresentantenString != "");
	
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
			if (IsSameType<TSequence, Peptide>::VALUE){
				GScorePeptideQGram(score, queryRecord, referenceRecord, scoringScheme, scoreMode);}
			else
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
	cout << "\tSearched "<<length(perfQueryAvgContainer)<<" queries and found "<<length(matchesCollection)<<" matches." << endl;
    PerformanceSample perfQueryAvg("Handle one query in GSearch average");
    //perfQueryAvg.takeAverageValues(perfQueryAvgContainer);
	//appendValue(performance, perfQueryAvg);
    perfSearchFunc.end();
	appendValue(performance, perfSearchFunc);
    return 0;
}

#endif  // GINGER_GSEARCH_GSEARCH_BASE_H_