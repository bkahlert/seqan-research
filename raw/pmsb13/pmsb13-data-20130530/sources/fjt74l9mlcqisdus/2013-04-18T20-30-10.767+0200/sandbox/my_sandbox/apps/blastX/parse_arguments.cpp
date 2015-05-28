#include "own_functions.h"

void dafault_values(Variable & comVal){
	comVal.fasta_file = "database.fa";
	comVal.fastq_file = "reads.fq"
	comVal.seed = 15; 
	comVal.size_alp = 5;
	comVal.numb_alp = 3;
}

int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal){
	// Setup ArgumentParser.
    ArgumentParser parser("blastX");

    addOption(parser, seqan::ArgParseOption("d", "DATABASE", "PATH OF THE TARGET DATABASE FASTA FILE",
        ArgParseArgument::STRING, "STRING"));

	addOption(parser, seqan::ArgParseOption("r", "READS", "PATH OF THE TARGET READ FASTQ FILE",
        ArgParseArgument::STRING, "STRING"));

	addOption(parser, seqan::ArgParseOption("s", "SEED", "SIZE OF THE PIECES FROM THE QGRAM INDEX",
        ArgParseArgument::INT, "INT"));

	addOption(parser, seqan::ArgParseOption("a", "SIZE_ALPHABET", "SIZE OF THE REDUCED AMINO ACID ALPHABET",
        ArgParseArgument::INT, "INT"));

	addOption(parser, seqan::ArgParseOption("n", "NUMBER_ALPHABET", "NUMBER OF THE DIFFERENT USED ALPHABETS",
        ArgParseArgument::INT, "INT"));

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
	getOptionValue(comVal.fastq_file, parser, "READS");
	getOptionValue(comVal.seed, parser, "SEED");
	getOptionValue(comVal.size_alp, parser, "SIZE_ALPHABET");
	getOptionValue(comVal.numb_alp, parser, "NUMBER_ALPHABET");
	return 0;
}