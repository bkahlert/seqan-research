// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

/**
*@brief parse_arguments.cpp
* Beinaltet Funktionen, die die Eingabeparameter ueber die Kommandozeile parst, checkt und evtl. default Werte einsetzt.
*/

#include "own_functions.h"


int CHECK_VALUES(Variable & comVal,StringSet<String<Dna> > & reads){
	/**@brief check_values Ueberprueft die Eingabeparameter auf ihre Richtigkeit
	*@param comVal Gespeicherte Eingabeparameter
	*@param reads Aus der Datei ausgelesene Reads  
	*@return Gibt 0 zurueck falls es erfolgreich war, ansonsten 1 wenn das parsen fehlschlaegt
	*/
	if (comVal.seed<2){
		comVal.seed = 2;
		cerr << "SEED too low put high to : "<< comVal.seed;
	}
	for (unsigned i = 0; i < length(reads); ++i){
		if (comVal.seed>(length(reads[i])/3) && length(reads[i])>=6){	//extra Abfrage:length(Reads[i])>=6) ->Programm abbrechen
			comVal.seed = (length(reads[i])/3);
			cerr << "SEED too high put down to : "<< comVal.seed;
		}
	}
	
	if (comVal.size_alp < 7 || comVal.size_alp > 20){
		cerr << "SIZE_ALPHABET MUST BE BETWEEN 7-20 " <<endl;
		return 1;
	}
	
	/*if (comVal.numb_alp<0){
		cerr << "NUMBER_ALPHABET MUST BE >0 "<<endl;
		return 1;
	}*/
		
	//if (comVal.hamming_distance<0 || comVal.hamming_distance>comVal.seed){
	if (comVal.hamming_distance>comVal.seed){
		cerr << "DISTANCE MUST BE BETWEEN 0 AND SEED "<<endl;
		return 1;
	}

	return 0;
}

void DEFAULT_VALUES(Variable & comVal){
	/**
	*@brief Initialisierung der Klasse Variable mit default Werten
	*@param comVal speichert default Werte
	*@return Kein Rueckgabewert
	*/
	comVal.fasta_file = "database.fa";
	comVal.fastq_file = "read.fq";
	comVal.seed = 15; 
	comVal.size_alp = 10;
	appendValue(comVal.numb_alp,49);
	comVal.hamming_distance = 0;
}


int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal){
	/**@brief Werte aus der Kommandozeile werden geparst 
	*@param argc Anzahl an Argumenten, die der Kommandozeile uebergeben werden
	*@param argv Speichert alle Argumente, die der Kommandozeile uebergeben werden
	*@param comVal Speichert die geparsten Werte 
	*@return Gibt 0 zurueck falls es erfolgreich war, ansonsten 1 wenn das parsen fehlschlaegt 
	*/
	
	//Optimierung der Kommandozeileneingabe
    ArgumentParser parser("blastX");

    addOption(parser, seqan::ArgParseOption("d", "DATABASE", "FILENAME OF THE TARGET DATABASE FASTA FILE",
        ArgParseArgument::STRING, "STRING"));

	addOption(parser, seqan::ArgParseOption("r", "READS", "FILENAME OF THE TARGET READ FASTQ FILE",
        ArgParseArgument::STRING, "STRING"));

	addOption(parser, seqan::ArgParseOption("s", "SEED", "SIZE OF A SHORT SUBSEQUENCE",
        ArgParseArgument::INTEGER, "INT"));

	addOption(parser, seqan::ArgParseOption("a", "SIZE_ALPHABET", "SIZE OF THE REDUCED AMINO ACID ALPHABET",
        ArgParseArgument::INTEGER, "INT"));

	addOption(parser, seqan::ArgParseOption("n", "NUMBER_ALPHABET", "NUMBER OF THE DIFFERENT USED ALPHABETSCORES: (1)Blosum30, (2)Blosum45, (3)Blosum80, (4)Pam120, (5)Pam200, (6)Pam250, (7)Pam40, (8)Vtml200",
		ArgParseArgument::STRING, "STRING"));
	
	addOption(parser, seqan::ArgParseOption("k", "DISTANCE", "NUMBER OF THE HAMMING DISTANCE BETWEEN PATTERN AND DATABASE",
		ArgParseArgument::INTEGER, "INT"));

    //Parse Kommandozeile
    ArgumentParser::ParseResult res = parse(parser, argc, argv); 

 
    //Abfrage ob das parsen erfolgreich ist 
	if (res != ArgumentParser::PARSE_OK){
		//ArgumentParser::PARSE_ERROR;
		return 1;
	}
        
    //Speicherung der Werte in Klasse comVal
	getOptionValue(comVal.fasta_file, parser, "DATABASE");
	getOptionValue(comVal.fastq_file, parser, "READS");
	getOptionValue(comVal.seed, parser, "SEED");
	getOptionValue(comVal.size_alp, parser, "SIZE_ALPHABET");
	getOptionValue(comVal.numb_alp, parser, "NUMBER_ALPHABET");
	getOptionValue(comVal.hamming_distance, parser, "DISTANCE");
	
	cout << "ARGUMENTS_ARE_PARSING"<<endl;
	return 0;
}

