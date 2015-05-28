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


template<typename TPath>
int readCommandlineParameters(TPath & clusteredDatabase, TPath & clusterRepresentanten, TPath & queries, TPath & outFile, string & scoreMode, int & argc, char const ** & argv) {
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("GSearch");

    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"Database\" \"ClusterMaster\" \"Queries\" \"outPathMatches\"");
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
			  "Mode of scoring function.\n\tLOCAL_ALIGNMENT_MODE\n\tCOMMON_QGRAM_MODE",
		seqan::ArgParseArgument::STRING, "STRING"));



    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

    getArgumentValue(clusteredDatabase, parser, 0);

    getArgumentValue(clusterRepresentanten, parser, 1);

    getArgumentValue(queries, parser, 2);

    getArgumentValue(outFile, parser, 3);
	
	getOptionValue(scoreMode, parser, "scoringMode");
    return 0;
}



int checkStreams(SequenceStream &myCluster, SequenceStream &myMaster, SequenceStream &myQueries, fstream &myOutput) {
    if (!isGood(myCluster)) {
        std::cerr << "ERROR: Could not open the clustered database file in function GSearch.\n";
        return 1;
    }

    if (!isGood(myMaster)) {
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
int GSearch(String<GMatch<TScore> > matches, String<PerformanceSample> & performance, int argc, char const ** argv) {
	PerformanceSample perfSearchFunc("GSearch Funktion");
    typedef Score<int> TScoringScheme;
	typedef Rna5String TSequence;

    String<char> queries;
    String<char> outFile;
    String<char> clusterRepresentanten;
	String<char> clusteredDatabase;
	string scoreMode;
    
	int res = readCommandlineParameters(clusteredDatabase, clusterRepresentanten, queries, outFile, scoreMode, argc, argv);
	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;
	
    SequenceStream myCluster(toCString(clusteredDatabase));
    SequenceStream myMaster(toCString(clusterRepresentanten));
    SequenceStream myQueries(toCString(queries));
    std::fstream myOutput(toCString(outFile), std::ios::binary | std::ios::out);
    if (checkStreams(myCluster, myMaster, myQueries, myOutput))
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
		
        while(!readRecord(referenceRecord.id, referenceRecord.seq, myMaster)) {
			GScore(score, queryRecord, referenceRecord, scoringScheme, scoreMode);
			GEvaluateScore(scoreStorage, score, referenceRecord);
		}
        open(myCluster, toCString(clusteredDatabase));
        while(!readRecord(clusterRecord.id, clusterRecord.seq, myCluster)) {
			if (GCheckClusterMaster(scoreStorage, clusterRecord)) {
				GComputeMatch(matches, queryRecord, clusterRecord, scoringScheme, scoreMode);
            }
        }
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