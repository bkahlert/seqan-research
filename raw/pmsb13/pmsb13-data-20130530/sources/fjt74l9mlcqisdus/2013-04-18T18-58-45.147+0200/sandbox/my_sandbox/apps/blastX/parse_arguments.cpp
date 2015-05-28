#include <seqan/arg_parse.h>
#include "own_functions.h"
using namespace std;
using namespace seqan;

struct Values{
	String<char> fasta_file;
	Values() : fasta_file("database.fa") {}
	String<char> fastq_file;
	Values() : fastq_file("reads.fq") {}
};

int PARSE_ARGUMENTS(int argc,char const ** argv,Values comVal){
	// Setup ArgumentParser.
    ArgumentParser parser("blastX");

    addOption(parser, seqan::ArgParseOption("d", "DATABASE", "PATH OF THE TARGET DATABASE FASTA-FILE",
        ArgParseArgument::STRING<char>, "STRING<char>"));

	 addOption(parser, seqan::ArgParseOption("r", "READS", "PATH OF THE TARGET READ FASTQ-FILE",
        ArgParseArgument::STRING<char>, "STRING<char>"));

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
	getOptionValue(comVal.fasta_file, parser, "READS");
	return 0;
}
