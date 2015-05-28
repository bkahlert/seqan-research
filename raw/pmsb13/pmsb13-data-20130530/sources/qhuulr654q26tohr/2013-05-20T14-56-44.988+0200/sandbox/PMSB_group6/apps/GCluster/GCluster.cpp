//Autor:Jakob

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

#include <seqan/GCluster.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	
	CharString inPath;
	CharString outPath;
	CharString outMasterPath;
	CharString sequenceFormat;
	int threshold;
	double lengthThreshold;
	CharString scoreMode;
	
	// Setup ArgumentParser.
	seqan::ArgumentParser parser("GCluster");
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
			 seqan::ArgParseOption("m", "scoringMode",
								   "Mode of scoring function.\n\tLOCAL_ALIGNMENT_MODE\n\tGLOBAL_ALIGNMENT_MODE\n\tCOMMON_QGRAM_MODE",
						  seqan::ArgParseArgument::STRING, "STRING"));
			addOption(
				parser,
			 seqan::ArgParseOption("l", "clusterLengthThreshold",
								   "Maximum relative length difference of the master sequence to the other sequences of the cluster.",
						  seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
			addOption(
				parser,
			 seqan::ArgParseOption("f", "sequenceFormat",
					  "Sequence format of database.\tDNA\tRNA\tAMINOACID",
			 seqan::ArgParseArgument::STRING, "STRING"));
			
			addDefaultValue(parser, "scoringMode", "COMMON_QGRAM_MODE");
			addDefaultValue(parser, "threshold", 50);
			addDefaultValue(parser, "clusterLengthThreshold", 0.5);
			addDefaultValue(parser, "sequenceFormat", "DNA");
			
			
			// Parse command line.
			seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
			
			// If parsing was not successful then exit with code 1 if there were errors.
			// Otherwise, exit with code 0 (e.g. help was printed).
			if (res != seqan::ArgumentParser::PARSE_OK)
				return res == seqan::ArgumentParser::PARSE_ERROR;
			
			getArgumentValue(inPath, parser, 0);
			getArgumentValue(outPath, parser, 1);
			getArgumentValue(outMasterPath, parser, 2);
			
			getOptionValue(threshold, parser, "threshold");
			getOptionValue(scoreMode, parser, "scoringMode");
			getOptionValue(lengthThreshold, parser, "clusterLengthThreshold");
			getOptionValue(sequenceFormat, parser, "sequenceFormat");
			
			string scoreModeString = toCString(scoreMode);
			
	String<PerformanceSample> performance;
	
	if (sequenceFormat=="DNA"){
		Dna5String sequenceFormatExample;
		int res = GCluster(performance, toCString(inPath), toCString(outPath), toCString(outMasterPath), sequenceFormatExample, threshold, lengthThreshold, scoreModeString);
		if (res!=0) {
			cerr << "Error in function GCluster" << endl;
			return 1;
		}
		return 0;
	}
	if (sequenceFormat=="RNA"){
		Rna5String sequenceFormatExample;
		int res = GCluster(performance, toCString(inPath), toCString(outPath), toCString(outMasterPath), sequenceFormatExample, threshold, lengthThreshold, scoreModeString);
		if (res!=0) {
			cerr << "Error in function GCluster" << endl;
			return 1;
		}
		return 0;
	}
	if (sequenceFormat=="AMINOACID"){
		Peptide sequenceFormatExample;
		int res = GCluster(performance, toCString(inPath), toCString(outPath), toCString(outMasterPath), sequenceFormatExample, threshold, lengthThreshold, scoreModeString);
		if (res!=0) {
			cerr << "Error in function GCluster" << endl;
			return 1;
		}
		return 0;
	}
	cerr << "Unknown sequence format: " << sequenceFormat << endl;
	return 1;
}