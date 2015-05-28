#include "own_functions.h"
#include <iostream>

using namespace std;
using namespace seqan;

struct Values{
	String<char> fasta_file;
	Values() : DNA_LENGTH("database.fa") {}
};

int PARSE_ARGUMENTS(int argc,char const ** argv,Values comVal){
	// Setup ArgumentParser.
    ArgumentParser parser("blastX");

    addOption(parser, seqan::ArgParseOption("d", "DATABASE", "PATH OF THE TARGET DATABASE FASTA FILE",
        ArgParseArgument::STRING, "STRING"));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
	if (res != ArgumentParser::PARSE_OK){
		ArgumentParser::PARSE_ERROR;
		return 1;
	}
        
    // Extract option values and print them.
	getOptionValue(comVal.fasta_file, parser, "DATABASE");
	return 0;
}