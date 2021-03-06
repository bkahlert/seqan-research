// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

#include "own_functions.h"

// INITIALISIERUNG DER KLASSE VARIABLE MIT
// DEFAULT WERTEN
void DEFAULT_VALUES(Variable & comVal){
	comVal.fasta_file = "database.fa";
	comVal.fastq_file = "read.fq";
	comVal.sam_file = "matches.sam";
	comVal.seed = 15; 
	comVal.size_alp = 5;
	comVal.numb_alp = 3;
}

// WERTE AUS KOMMANDOZEILE WERDEN GEPARST UND IN comVal GESCHRIEBEN
// WENN DAS PARSEN FEHL SCHLAEGT DANN RETURN WERT 1 
// SONST SCHLIE�t DIE FUNKTION MIT RETURN WERT 0
int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal){
	// ArgumentParser.
    ArgumentParser parser("blastX");

    addOption(parser, seqan::ArgParseOption("d", "DATABASE", "FILENAME OF THE TARGET DATABASE FASTA FILE",
        ArgParseArgument::STRING, "STRING"));

	addOption(parser, seqan::ArgParseOption("r", "READS", "FILENAME OF THE TARGET READ FASTQ FILE",
        ArgParseArgument::STRING, "STRING"));

	addOption(parser, seqan::ArgParseOption("o", "OUTPUT", "FILENAME OF THE OUTPUT SAM FILE",
        ArgParseArgument::STRING, "STRING"));

	addOption(parser, seqan::ArgParseOption("s", "SEED", "SIZE OF A SHORT SUBSEQUENCE",
        ArgParseArgument::INTEGER, "INT"));

	addOption(parser, seqan::ArgParseOption("a", "SIZE_ALPHABET", "SIZE OF THE REDUCED AMINO ACID ALPHABET",
        ArgParseArgument::INTEGER, "INT"));

	addOption(parser, seqan::ArgParseOption("n", "NUMBER_ALPHABET", "NUMBER OF THE DIFFERENT USED ALPHABETS",
		ArgParseArgument::INTEGER, "INT"));

	

    // Parse Commandozeile
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // wenn das Parsen fehl schlaegt dann return - Wert 1 
    // sonst schlie�t die Funktion mit return - Wert 0
	if (res != ArgumentParser::PARSE_OK){
		ArgumentParser::PARSE_ERROR;
		return 1;
	}
        
    // Speicherung der Werte in Klasse comVal
	getOptionValue(comVal.fasta_file, parser, "DATABASE");
	getOptionValue(comVal.fastq_file, parser, "READS");
	getOptionValue(comVal.sam_file, parser, "OUTPUT");
	getOptionValue(comVal.seed, parser, "SEED");
	getOptionValue(comVal.size_alp, parser, "SIZE_ALPHABET");
	getOptionValue(comVal.numb_alp, parser, "NUMBER_ALPHABET");
	
	return 0;
}

