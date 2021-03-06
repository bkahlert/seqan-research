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
	setShortDescription(parser, "Searching a list of queries in a (clustered) database.");
	setVersion(parser, "1.0");
	setDate(parser, "May 2013");
	addDescription(parser,"Searching a query in a large database can be accelerated if the sequences are first assign to a number of clusters representing similar sequences. The query is then first matched against the cluster masters and then against every sequence in the cluster. GSearch reads a fasta file containing the queries and searchs each query in a database. If a file containing master sequences is specified it assumes that the database is clustered and checks only the sequences of the cluster whose master scored highest.");

    addUsageLine(parser,
				 "[\\fIOPTIONS\\fP] \"Database_InPath\" \"Queries_InPath\" \"Matches_OutPath\"");
    addArgument(
        parser,
		seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE,
                                "Path to input cluster FASTA"));
    addArgument(
        parser,
		seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE,
                                "Path to input queries FASTA"));
    addArgument(
        parser,
		seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUTFILE,
                                "Path to output file FASTA"));
    addSection(parser, "Modification Options");
    addOption(
        parser,
		seqan::ArgParseOption("M", "masterSequences","Path to fasta file containing master sequences",seqan::ArgParseArgument::INPUTFILE,"STRING"));
    addOption(
        parser,
        seqan::ArgParseOption("m", "scoringMode",
                              "Mode of scoring function.",
                              seqan::ArgParseArgument::STRING, "STRING"));
    addOption(
        parser,
        seqan::ArgParseOption("t", "threshold",
                              "Minimum score of reported match",
							  seqan::ArgParseArgument::INTEGER, "INTEGER"));
	addOption(
		parser,
		seqan::ArgParseOption("f", "sequenceFormat",
							  "Sequence format of database.",
						seqan::ArgParseArgument::STRING, "STRING"));

    addDefaultValue(parser, "scoringMode", "COMMON_QGRAM_MODE");
    addDefaultValue(parser, "masterSequences", "NOFILE");
	addDefaultValue(parser, "threshold", "1");
	addDefaultValue(parser, "sequenceFormat", "DNA");
	setMinValue(parser, "threshold", "0");
	setMaxValue(parser, "threshold", "1000");
	setValidValues(parser, "scoringMode", "COMMON_QGRAM_MODE LOCAL_ALIGNMENT_MODE GLOBAL_ALIGNMENT_MODE");
	setValidValues(parser, "sequenceFormat", "DNA RNA AMINOACID");
	
	addText(parser,"LOCAL_ALIGNMENT_MODE suming up local alignments.");
	addText(parser,"GLOBAL_ALIGNMENT_MODE computes a global alignment with no gap costs.");
	addText(parser,"COMMON_QGRAM_MODE counts common q-grams (DNA/RNA=8;AMINOACID=3). The score is the fraction of q-grams of the query found in the target multiplied with 1000. [recommended mode]");
	addTextSection(parser, "Examples");
	addText(parser,"The mentioned test files can be found in \\fBrepository\\fP/tests/autoTestData/unmutable/");
	addListItem(parser,"GSearch \\fBclustered.pfam.CQ40.fasta queries_pfam.included.fasta Matches_OutPath\\fP -M \\fBmaster.pfam.CQ40.fasta\\fP -f AMINOACID -t 200", "Matches have at least 20% similar 3-grams to the target. Uses a selection of the pfam data where the three queries are included in the database.");
	addListItem(parser,"GSearch \\fBclustered.rfam.CQ20.fasta queries_rfam.included.fasta Matches_OutPath\\fP -M \\fBmaster.rfam.CQ20.fasta\\fP -f RNA -t 200", "Matches have at least 20% similar 8-grams to the cluster master. Uses a selection of the rfam data where the three queries are included in the database.");
	


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