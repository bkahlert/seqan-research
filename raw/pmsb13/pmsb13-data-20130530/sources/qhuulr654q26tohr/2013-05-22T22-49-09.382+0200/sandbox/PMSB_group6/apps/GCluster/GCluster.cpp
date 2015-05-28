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
	setShortDescription(parser, "Clustering of a fasta database.");
	setVersion(parser, "1.0");
	setDate(parser, "May 2013");
	addDescription(parser,"Searching a query in a large database can be accelerated if the sequences are first assign to a number of clusters representing similar sequences. The query is then first matched against the cluster masters and then against every sequence in the cluster. GCluster reads a database in fasta format and creates a new clustered database as well as a file storing all master sequences.");
    addUsageLine(parser,"[\\fIOPTIONS\\fP] \"Database_InPath\" \"ClusteredDatabase_OutPath\" \"ClusterMasters_OutPath\"");
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE,
                                "Path to input fasta"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUTFILE,
                                "Path to database output fasta"));
    addArgument(
        parser,
        seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUTFILE,
                                "Path to master sequences output fasta"));
    addSection(parser, "Modification Options");
    addOption(
        parser,
        seqan::ArgParseOption("t", "threshold",
                              "threshold for clustering local alignment score.",
                              seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(
        parser,
        seqan::ArgParseOption("m", "scoringMode",
                              "Mode of scoring function.",
							  seqan::ArgParseArgument::STRING, "STRING"));
    addOption(
        parser,
        seqan::ArgParseOption("l", "clusterLengthThreshold",
                              "Maximum relative length difference of the master sequence to the other sequences of the cluster.",
                              seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(
        parser,
        seqan::ArgParseOption("f", "sequenceFormat",
                              "Sequence format of database.",
                              seqan::ArgParseArgument::STRING, "STRING"));

    addDefaultValue(parser, "scoringMode", "COMMON_QGRAM_MODE");
    addDefaultValue(parser, "threshold", "50");
    addDefaultValue(parser, "clusterLengthThreshold", "0.5");
    addDefaultValue(parser, "sequenceFormat", "DNA");
    setMinValue(parser, "threshold", "0");
    setMaxValue(parser, "threshold", "1000");
    setMinValue(parser, "clusterLengthThreshold", "0");
    setMaxValue(parser, "clusterLengthThreshold", "1");
    setValidValues(parser, "scoringMode", "COMMON_QGRAM_MODE LOCAL_ALIGNMENT_MODE GLOBAL_ALIGNMENT_MODE");
    setValidValues(parser, "sequenceFormat", "DNA RNA AMINOACID");
	
	addText(parser,"LOCAL_ALIGNMENT_MODE suming up local alignments.");
	addText(parser,"GLOBAL_ALIGNMENT_MODE computes a global alignment with no gap costs.");
	addText(parser,"COMMON_QGRAM_MODE counts common q-grams (DNA/RNA=8;AMINOACID=3). The score is the fraction of q-grams of the query found in the target multiplied with 1000. [recommended]");
	addTextSection(parser, "Examples");
	addText(parser,"The mentioned test files can be found in \\fBrepository\\fP/tests/autoTestData/unmutable/");
	addListItem(parser,"GCluster \\fBpfam.100K.fasta ClusteredDatabase_OutPath ClusterMasters_OutPath\\fP -f AMINOACID -t 80", "Cluster members have 8% similar 3-grams to the cluster master. Uses a selection of the pfam data.");
	addListItem(parser,"GCluster \\fBrfam.10K.fasta ClusteredDatabase_OutPath ClusterMasters_OutPath\\fP -f RNA -t 80 -l 0.1", "Cluster members have 8% similar 8-grams and their length varies only by 10% from the cluster master. Uses a selection of the rfam data.");

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

    if (sequenceFormat=="DNA") {
        Dna5String sequenceFormatExample;
        int res = GCluster(performance, toCString(inPath), toCString(outPath), toCString(outMasterPath), sequenceFormatExample, threshold, lengthThreshold, scoreModeString);
        if (res!=0) {
            cerr << "Error in function GCluster" << endl;
            return 1;
        }
        return 0;
    }
    if (sequenceFormat=="RNA") {
        Rna5String sequenceFormatExample;
        int res = GCluster(performance, toCString(inPath), toCString(outPath), toCString(outMasterPath), sequenceFormatExample, threshold, lengthThreshold, scoreModeString);
        if (res!=0) {
            cerr << "Error in function GCluster" << endl;
            return 1;
        }
        return 0;
    }
    if (sequenceFormat=="AMINOACID") {
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