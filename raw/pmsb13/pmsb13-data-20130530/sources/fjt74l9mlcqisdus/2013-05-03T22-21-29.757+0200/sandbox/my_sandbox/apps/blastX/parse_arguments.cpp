// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

/**
*@brief Datei parse_arguments.cpp
* Beinaltet Funktionen, die die Eingabeparameter ueber die Kommandozeile parst.
*/

#include "own_functions.h"


int check_values(Variable & comVal,StringSet<String<Dna5>> & Reads){
	if (comVal.seed<2){
		comVal.seed = 2;
		cerr << "Seed_size too low put high to : "<< comVal.seed;
	}
	for (int i=0; i<length(Reads);++i){
		cout << length(Reads[i])<<"\t"<<(length(Reads[i])/3)<<endl;
		if (comVal.seed>(length(Reads[i])/3) && length(Reads[i])>=6){
			comVal.seed = (length(Reads[i])/3);
			cerr << "Seed_size too high put down to : "<< comVal.seed;
		}
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
	comVal.size_alp = 5;
	comVal.numb_alp = 3;
	comVal.hamming_distance = 0;
}


int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal){
	/**@brief Werte aus der Kommandozeile werden geparst
	*@param comVal Speichert die geparsten Werte  
	*@param argc Anzahl an Argumenten, die der Kommandozeile uebergeben werden
	*@param argv Speichert alle Argumente, die der Kommandozeile uebergeben werden
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

	addOption(parser, seqan::ArgParseOption("n", "NUMBER_ALPHABET", "NUMBER OF THE DIFFERENT USED ALPHABETS",
		ArgParseArgument::INTEGER, "INT"));
	
	addOption(parser, seqan::ArgParseOption("k", "DISTANCE", "NUMBER OF THE HAMMING DISTANCE BETWEEN PATTERN AND DATABASE",
		ArgParseArgument::INTEGER, "INT"));

    //Parse Kommandozeile
    ArgumentParser::ParseResult res = parse(parser, argc, argv); 

 
    //Abfrage ob das parsen erfolgreich ist 
	if (res != ArgumentParser::PARSE_OK){
		ArgumentParser::PARSE_ERROR;
		return 1;
	}
        
    //Speicherung der Werte in Klasse comVal
	getOptionValue(comVal.fasta_file, parser, "DATABASE");
	getOptionValue(comVal.fastq_file, parser, "READS");
	getOptionValue(comVal.seed, parser, "SEED");
	getOptionValue(comVal.size_alp, parser, "SIZE_ALPHABET");
	getOptionValue(comVal.numb_alp, parser, "NUMBER_ALPHABET");
	getOptionValue(comVal.hamming_distance, parser, "DISTANCE");
	
	cerr << "ARGUMENTS_ARE_PARSING"<<endl;
	return 0;
}

