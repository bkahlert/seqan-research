#ifndef GINGER_SEARCH_SEARCH_BASE_H_
#define GINGER_SEARCH_SEARCH_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <fstream>
#include <seqan/structs/GFastaRecord.h>
#include <seqan/structs/GMatch.h>
#include <seqan/structs/GScoreStorage.h>
#include <seqan/structs/PerformanceSample.h>
#include <seqan/structs/MemorySample.h>
#include <seqan/Search/GScore.h>
#include <seqan/Search/GEvaluateScore.h>
#include <seqan/Search/GParseMasterId.h>
#include <seqan/Search/GCheckClusterMaster.h>
#include <seqan/Search/GComputeMatch.h>
#include <seqan/Search/GPostProcessMatches.h>
#include <seqan/Search/GWriteMatches.h>

using namespace seqan;
using namespace std;


template<typename TPath>
int readCommandlineParameters(TPath & clusteredDatabase, TPath & clusterRepresentanten, TPath & queries, TPath & outFile, int & argc, char const ** & argv) {
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("search");

    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"Database\" \"ClusterMaster\" \"Queries\" \"outPathMatches\"");
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to input cluster FASTA"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to input Master FASTA"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to input queries FASTA"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to output file FASTA"));



    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    getArgumentValue(clusteredDatabase, parser, 0);

    getArgumentValue(clusterRepresentanten, parser, 1);

    getArgumentValue(queries, parser, 2);

    getArgumentValue(outFile, parser, 3);
    return 0;
}



int checkStreams(SequenceStream &myCluster, SequenceStream &myMaster, SequenceStream &myQueries, fstream &myOutput) {
    if (!isGood(myCluster)) {
        std::cerr << "ERROR: Could not open the Cluster file.\n";
        return 1;
    }

    if (!isGood(myMaster)) {
        std::cerr << "ERROR: Could not open the MasterSequences file.\n";
        return 1;
    }
    if (!isGood(myQueries)) {
        std::cerr << "ERROR: Could not open the Queries file.\n";
        return 1;
    }
    if (!myOutput.good()) {
        std::cerr << "ERROR: Could not open the output file.\n";
        return 1;
    }
    return 0;
}



int Search(String<PerformanceSample> & performance, int argc, char const ** argv) {
	PerformanceSample perfSearchFunc("Search Funktion");
    typedef Score<int> TScoringScheme;
    typedef int TScore;
	typedef DnaString TSequence;

    String<char> queries;
    String<char> outFile;
    String<char> clusterRepresentanten;
    String<char> clusteredDatabase;
    if (readCommandlineParameters(clusteredDatabase, clusterRepresentanten, queries, outFile, argc, argv))
        return 1;
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
	String<GMatch<TScore> > matches;
	String<PerformanceSample> perfQueryAvgContainer;

	while (!readRecord(queryRecord.id, queryRecord.seq, myQueries)) {
		PerformanceSample perfQueryAvg("");
		open(myMaster, toCString(clusterRepresentanten));
		GScoreStorage<TScore> scoreStorage;
		
        while(!readRecord(referenceRecord.id, referenceRecord.seq, myMaster)) {
			GScore(score, queryRecord, referenceRecord); //, scoringScheme
			GEvaluateScore(scoreStorage, score, referenceRecord);
		}
        open(myCluster, toCString(clusteredDatabase));
        while(!readRecord(clusterRecord.id, clusterRecord.seq, myCluster)) {
			if (GCheckClusterMaster(scoreStorage, clusterRecord)) {
				GComputeMatch(matches, queryRecord, clusterRecord, scoringScheme);
            }
        }
        GPostProcessMatches(matches);
		GWriteMatches(myOutput, matches, queryRecord);
		
        perfQueryAvg.end();
		appendValue(perfQueryAvgContainer, perfQueryAvg);
    }
    PerformanceSample perfQueryAvg("Handle one query in Search average");
    perfQueryAvg.takeAverageValues(perfQueryAvgContainer);
	appendValue(performance, perfQueryAvg);
    perfSearchFunc.end();
	appendValue(performance, perfSearchFunc);
    return 0;
}

#endif  // GINGER_SEARCH_SEARCH_BASE_H_