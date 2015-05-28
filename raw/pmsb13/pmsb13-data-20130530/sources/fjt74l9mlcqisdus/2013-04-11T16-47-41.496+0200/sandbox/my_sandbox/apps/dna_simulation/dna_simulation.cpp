#include <iostream>
#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;


void PARSE_ARGUMENTS(int argc,char const ** argv){
	// Setup ArgumentParser.
    ArgumentParser parser("dna_simulation");

    addArgument(parser, seqan::ArgParseArgument(ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("l", "DNA_LENGTH", "Length of the building DNA.",
        ArgParseArgument::INTEGER, "INT"));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
        
    // Extract option values and print them.
    unsigned DNA_LENGTH = 0;
    getOptionValue(length, parser, "DNA_LENGTH");

    cout << "DNA_LENGTH\t" << period << endl;
}



int main(int argc, char const ** argv)
{
	// take the comandline arguments and save them in ...
	PARSE_ARGUMENTS(argc,argv);
    return 0;
}