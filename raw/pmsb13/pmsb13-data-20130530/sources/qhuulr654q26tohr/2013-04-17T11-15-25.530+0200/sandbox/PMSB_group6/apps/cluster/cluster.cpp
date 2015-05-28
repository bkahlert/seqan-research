#include <iostream>
#include <seqan/arg_parse.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("cluster");
    
    addArgument(parser, seqan::ArgParseArgument(
	seqan::ArgParseArgument::STRING, "Path to input fasta"));
    addArgument(parser, seqan::ArgParseArgument(
	seqan::ArgParseArgument::STRING, "Path to output fasta"));
    
    addOption(parser, seqan::ArgParseOption(
	"t", "threshold", "threshold for clustering local alignment score",
	seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
    
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
	return res == seqan::ArgumentParser::PARSE_ERROR;
    
    // Extract option values and print them.
	String<Char> inPath;
	String<Char> outPath;
	double threshold;
	getArgumentValue(inPath, parser, 0);
	getArgumentValue(outPath, parser, 1);
	getOptionValue(threshold, parser, "threshold");
	
	std::cout << "inPath   \t" << inPath << '\n'
	<< "outPath\t" << outPath << '\n'
	<< "threshold     \t" << threshold << '\n';
	
	return 0;
}