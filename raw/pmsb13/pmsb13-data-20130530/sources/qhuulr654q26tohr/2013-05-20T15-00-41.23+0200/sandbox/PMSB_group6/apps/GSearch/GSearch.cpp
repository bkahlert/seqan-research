//autor: Jakob

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

#include <seqan/GSearch.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {

	CharString queries;
	CharString outFile;
	CharString clusterRepresentanten;
	CharString clusteredDatabase;
	CharString sequenceFormat;
	CharString scoreMode;
	int threshold;

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
	addOption(
		parser,
		seqan::ArgParseOption("f", "sequenceFormat",
							  "Sequence format of database.\tDNA\tRNA\tAMINOACID",
						seqan::ArgParseArgument::STRING, "STRING"));

    addDefaultValue(parser, "scoringMode", "COMMON_QGRAM_MODE");
    addDefaultValue(parser, "masterSequences", "");
	addDefaultValue(parser, "threshold", 1);
	addDefaultValue(parser, "sequenceFormat", "DNA");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

	getArgumentValue(clusteredDatabase, parser, 0);
	getArgumentValue(queries, parser, 1);
	getArgumentValue(outFile, parser, 2);
	
	getOptionValue(scoreMode, parser, "scoringMode");
	getOptionValue(clusterRepresentanten, parser, "masterSequences");
	getOptionValue(threshold, parser, "threshold");
	getOptionValue(sequenceFormat, parser, "sequenceFormat");

	string scoreModeString = toCString(scoreMode);

	String<PerformanceSample> performance;
	int databaseLength;

    if (sequenceFormat=="DNA") {
        String<GMatch<int, Dna5String> > matches;
		int res = GSearch(matches, databaseLength, performance, toCString(clusteredDatabase), toCString(clusterRepresentanten), toCString(queries), toCString(outFile), threshold, scoreModeString);
        if (res!=0) {
			cerr << "Error in function GSearch" << endl;
            return 1;
        }
        return 0;
    }
    if (sequenceFormat=="RNA") {
		String<GMatch<int, Rna5String> > matches;
		int res = GSearch(matches, databaseLength, performance, toCString(clusteredDatabase), toCString(clusterRepresentanten), toCString(queries), toCString(outFile), threshold, scoreModeString);
        if (res!=0) {
			cerr << "Error in function GSearch" << endl;
            return 1;
        }
        return 0;
    }
    if (sequenceFormat=="AMINOACID") {
		String<GMatch<int, Peptide> > matches;
		int res = GSearch(matches, databaseLength, performance, toCString(clusteredDatabase), toCString(clusterRepresentanten), toCString(queries), toCString(outFile), threshold, scoreModeString);
        if (res!=0) {
			cerr << "Error in function GSearch" << endl;
            return 1;
        }
        return 0;
    }
    cerr << "Unknown sequence format: " << sequenceFormat << endl;
    return 1;
}