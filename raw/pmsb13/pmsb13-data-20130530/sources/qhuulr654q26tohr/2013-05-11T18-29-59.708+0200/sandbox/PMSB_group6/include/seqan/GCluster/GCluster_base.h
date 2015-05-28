#ifndef GINGER_CLUSTER_CLUSTER_BASE_H_
#define GINGER_CLUSTER_CLUSTER_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/PerformanceSample.h>
#include <seqan/GStructs/MemorySample.h>
#include <seqan/GSearch/GScore.h>
#include <seqan/GCluster/GAssignCluster.h>


using namespace seqan;
using namespace std;

template<typename TPath, typename TThres>
int readCommandlineParameters(TPath & inPath, TPath & outPath, TPath & outMasterPath,
							  string & scoreMode, TThres & threshold, int & argc, char const ** & argv) {
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("cluster");
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"Database\" \"outPathDatabase\" \"outPathClusterMasters\"");
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to input fasta"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to database output fasta"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
                                "Path to master sequences output fasta"));
    addSection(parser, "Modification Options");
	addOption(
		parser,
		seqan::ArgParseOption("t", "threshold",
							  "threshold for clustering local alignment score",
						seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(
		parser,
		seqan::ArgParseOption("m", "scoring mode",
							  "mode of scoring function.\n\tLOCAL_ALIGNMENT_MODE\n\tCOMMON_QGRAM_MODE",
						seqan::ArgParseArgument::STRING, "STRING"));
	
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    getArgumentValue(inPath, parser, 0);
	getArgumentValue(outPath, parser, 1);
	getArgumentValue(outMasterPath, parser, 2);
	getArgumentValue(scoreMode, parser, 3);
	
	getOptionValue(threshold, parser, "threshold");
	getOptionValue(scoreMode, parser, "scoring mode");
    return 0;
}

int checkStreams(SequenceStream &inStream, SequenceStream &outStream, SequenceStream &outMasterStream) {
    if (!isGood(inStream)) {
        std::cerr << "ERROR: Could not open the input file.\n";
        return 1;
    }
    if (!isGood(outStream)) {
        std::cerr << "ERROR: Could not open the output file.\n";
        return 1;
    }
    if (!isGood(outMasterStream)) {
        std::cerr << "ERROR: Could not open the output file.\n";
        return 1;
    }
    return 0;
}


int GCluster(String<PerformanceSample> & performance, int argc, char const ** argv) {
    PerformanceSample perfCluster("GCluster Total");
    typedef int TThreshold;
    typedef Score<int> TScoringScheme;
    typedef GFastaRecord<DnaString> TRecord;
    typedef Iterator<String<TRecord> , Rooted>::Type TMasterIt;
    typedef int TScore;

    String<char> inPath;
    String<char> outPath;
    String<char> outMasterPath;
	TThreshold threshold;
	string scoreMode;
    if (readCommandlineParameters(inPath, outPath, outMasterPath, threshold, scoreMode, argc, argv))
        return 1;
    SequenceStream inStream(toCString(inPath));
    SequenceStream outStream(toCString(outPath), SequenceStream::WRITE);
    SequenceStream outMasterStream(toCString(outMasterPath), SequenceStream::WRITE);
    if (checkStreams(inStream, outStream, outMasterStream))
        return 1;

    String<TRecord> masterSequences;
    TScoringScheme scoringScheme(1, -2, -1);
    TMasterIt masterIt;
    TRecord record;

    String<PerformanceSample> perfEntryAvgContainer;
    while (!readRecord(record.id, record.seq, inStream)) {
        PerformanceSample perfEntryAvg("");
        bool createNewCluster = true;
        masterIt = begin(masterSequences);
        for (; !atEnd(masterIt); goNext(masterIt)) {
            TScore score;
			GScore(score, record, *masterIt, scoringScheme, scoreMode);
            if (threshold <= score) {
                createNewCluster = false;
                assignToCluster(outStream, record, *masterIt);
                break;
            }
        }
        if (createNewCluster) {
            appendValue(masterSequences, TRecord(record.id, record.seq));
            assignToCluster(outStream, record, record);
        }
        perfEntryAvg.end();
        appendValue(perfEntryAvgContainer, perfEntryAvg);
    }
    for (masterIt=begin(masterSequences); !atEnd(masterIt); goNext(masterIt)) {
        writeRecord(outMasterStream, (*masterIt).id, (*masterIt).seq);
    }

    //MemorySample ms("End of cluster function");
    //ms.printToStdout();

    PerformanceSample perfEntryAvg("Handle one database entry in Cluster average");
    perfEntryAvg.takeAverageValues(perfEntryAvgContainer);
    appendValue(performance, perfEntryAvg);
    perfCluster.end();
    appendValue(performance, perfCluster);
    return 0;
}

#endif  // GINGER_CLUSTER_CLUSTER_BASE_H_