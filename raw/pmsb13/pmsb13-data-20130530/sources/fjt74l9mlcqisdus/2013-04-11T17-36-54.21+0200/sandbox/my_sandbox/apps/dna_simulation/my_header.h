#ifndef SANDBOX_MY_SANDBOX_APPS_DNA_SIMULATION_MY_HEADER_

#include <seqan/arg_parse.h>
using namespace std;
using namespace seqan;

struct Values{
	int DNA_LENGTH;
	Values() : DNA_LENGTH(100) {}
};



int PARSE_ARGUMENTS(int argc,char const ** argv,Values comVal){
	// Setup ArgumentParser.
    ArgumentParser parser("dna_simulation");

    addOption(parser, seqan::ArgParseOption("l", "DNA_LENGTH", "Length of the building DNA.",
        ArgParseArgument::INTEGER, "INT"));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
	if (res != ArgumentParser::PARSE_OK){
		ArgumentParser::PARSE_ERROR;
		return 1;
	}
        
    // Extract option values and print them.
	comVal.DNA_LENGTH = 42;
    getOptionValue(comVal.DNA_LENGTH, parser, "DNA_LENGTH");
	return 0;
}


#define SANDBOX_MY_SANDBOX_APPS_DNA_SIMULATION_MY_HEADER_


#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_DNA_SIMULATION_MY_HEADER_
