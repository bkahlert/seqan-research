#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_OURSEARCH_BASE_H_
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_OURSEARCH_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <seqan/structs/PerformanceSample.h>
#include <seqan/structs/ScoredSequence.h>
#include <seqan/structs/FastaRecord.h>
#include <seqan/Search/getMasterId.h>
#include <seqan/Search/GScore.h>

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

bool myVergleich (ScoredSequence i,ScoredSequence j) {
    return (i.score<j.score);
}

int sortiereMatches (String<ScoredSequence> &unsortetResults) {

    std::sort (begin(unsortetResults), end(unsortetResults), myVergleich);
    return 0;
}



int Search(String<PerformanceSample> & performance, int argc, char const ** argv) {
	PerformanceSample perfSearchFunc("Search Funktion");
    typedef int TThreshold;
    typedef Score<int> TScoringScheme;
    typedef FastaRecord<RnaString> TRecord;
    typedef Iterator<String<TRecord> , Rooted>::Type TMasterIt;
    typedef int TScore;

    String<char> masterId;
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

    String<TRecord> masterSequences;
    TScoringScheme scoringScheme(1, -2, -1);
    TMasterIt masterIt;
    TRecord queryRecord;
    TRecord referenceRecord;
    TRecord clusterRecord;

    String<char> maxId;
    int score;
    ScoredSequence PartResults;
    String<ScoredSequence> Results;

	String<PerformanceSample> perfQueryAvgContainer;
	while (!readRecord(queryRecord.id, queryRecord.seq, myQueries)) {
		PerformanceSample perfQueryAvg("");
        open(myMaster, toCString(clusterRepresentanten));
        int maxScore=-1000;//nur provisorisch... muss überarbeitet werden
        clear(Results);
		
        while(!readRecord(referenceRecord.id, referenceRecord.seq, myMaster)) {
            GScore(score, queryRecord, referenceRecord, scoringScheme);
			// GEvaluateScore function
			// GScoreStorage struct
            if (maxScore <= score) {
                maxScore=score;
                maxId=referenceRecord.id;
            }
        }
        open(myCluster, toCString(clusteredDatabase));
        while(!readRecord(clusterRecord.id, clusterRecord.seq, myCluster)) {
			// if GCheckClusterMaster( u.a. GScoreStorage)
			//		GComputeMatch(String<GMatch>, ...)
			// GMatch struct
            masterId=getMasterId(clusterRecord.id);
            if (masterId==maxId) {
                GScore(score, queryRecord, clusterRecord, scoringScheme);
                PartResults.score = score;
                PartResults.id=clusterRecord.id;
                appendValue(Results, PartResults);
            }
        }

		// GPostProcessMatches(String<GMatch>)
        sortiereMatches(Results);
        stringstream ss;

		// GWriteMatches(...)
        seqan::CharString buffer;
        append(buffer, "\n Query ID =");
        append(buffer, queryRecord.id);
        for(int i=0; i<length(Results); i++) {
            append(buffer, "\n \t MasterSequence ID =");
            append(buffer, Results[i].id);
            append(buffer, "\t AlignmentScore =");
            ss.str("");
            ss.clear();
            ss << Results[i].score;
            append(buffer, ss.str());
        }

        if(seqan::streamWriteBlock(myOutput, &buffer[0], length(buffer)) != length(buffer)) {
            std::cerr << "ERROR: Could not print output.\n";
            return 1;
        }
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

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_OURSEARCH_BASE_H_