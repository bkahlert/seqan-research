#ifndef GINGER_GSEARCH_GSEARCH_BASE_H_
#define GINGER_GSEARCH_GSEARCH_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <fstream>
#include "GStructs/GFastaRecord.h"
#include "GStructs/GMatch.h"
#include "GStructs/GScoreStorage.h"
#include "GStructs/PerformanceSample.h"
#include "GStructs/MemorySample.h"
#include "GSearch/GScore.h"
#include "GSearch/GEvaluateScore.h"
#include "GSearch/GParseMasterId.h"
#include "GSearch/GCheckClusterMaster.h"
#include "GSearch/GComputeMatch.h"
#include "GSearch/GPostProcessMatches.h"
#include "GSearch/GWriteMatches.h"

using namespace seqan;
using namespace std;


template<typename TPath, typename TScore>
int readCommandlineParameters(TPath & clusteredDatabase, TPath & clusterRepresentanten, TPath & queries, TPath & outFile, string & scoreMode, TScore & threshold, int & argc, char const ** & argv) {
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("GSearch");

    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"Database\" \"Queries\" \"outPathMatches\"");
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to input cluster FASTA"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to input queries FASTA"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
								"Path to output file FASTA"));
	addSection(parser, "Modification Options");
	addOption(
		parser,
		seqan::ArgParseOption("M", "masterSequences","Path to input Master FASTA",seqan::ArgParseArgument::STRING,"STRING"));
	addOption(
		parser,
		seqan::ArgParseOption("m", "scoringMode",
							  "Mode of scoring function.\tLOCAL_ALIGNMENT_MODE\tCOMMON_QGRAM_MODE",
						seqan::ArgParseArgument::STRING, "STRING"));
	addOption(
		parser,
		seqan::ArgParseOption("t", "threshold",
							  "Minimum score of reported Match",
						seqan::ArgParseArgument::INTEGER, "INTEGER"));
	
	addDefaultValue(parser, "scoringMode", "COMMON_QGRAM_MODE");
	addDefaultValue(parser, "masterSequences", "");
	addDefaultValue(parser, "threshold", 1);

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

    getArgumentValue(clusteredDatabase, parser, 0);

    getArgumentValue(queries, parser, 1);

    getArgumentValue(outFile, parser, 2);
	
	getOptionValue(scoreMode, parser, "scoringMode");
	getOptionValue(clusterRepresentanten, parser, "masterSequences");
	getOptionValue(threshold, parser, "threshold");
    return 0;
}



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


template <typename TScore>
int GSearch(String<GMatch<TScore> > & matchesCollection, int & dataBaseLength, String<PerformanceSample> & performance, int argc, char const ** argv) {
	PerformanceSample perfSearchFunc("GSearch Funktion");
    typedef Score<int> TScoringScheme;
	typedef Rna5String TSequence;

    String<char> queries;
    String<char> outFile;
    String<char> clusterRepresentanten;
	String<char> clusteredDatabase;
	string scoreMode;
	TScore threshold;
    
	int res = readCommandlineParameters(clusteredDatabase, clusterRepresentanten, queries, outFile, scoreMode, threshold, argc, argv);
	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;
	
	bool MASTER_FILE_GIVEN = (clusterRepresentanten != "");
	
    SequenceStream myCluster(toCString(clusteredDatabase));
    SequenceStream myMaster(toCString(clusterRepresentanten));
    SequenceStream myQueries(toCString(queries));
    std::fstream myOutput(toCString(outFile), std::ios::binary | std::ios::out);
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
		String<GMatch<TScore> > matches;
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