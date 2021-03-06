// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

#include "own_functions.h"

// MAIN FUNKTION -------------------------------------------------------------------------------
int main(int argc, char const ** argv){
	/*
	// Klasse Variable enth�lt alle Parameter die der
	// Benutzer �ber Kommandozeile festlegen kann
	// -> own_functions.h
	Variable comVal;
	// comVal wird mit default Parameter initialisiert
	// -> parse_arguments.cpp
	DEFAULT_VALUES(comVal);
	// Eingabe vom Benutzer wird gelesen und default
	// Parameter wenn Eingabe vorhanden ist ueberschrieben
	if (argc>=2){
		// -> parse_arguments.cpp
		if (PARSE_ARGUMENTS(argc,argv,comVal))
			return 1;
	}
	// Alphabet enth�lt nach GET_ALPHABET numb_alp mal 
	// String<int> der laenge einundzwanzig wobei jede Position
	// eine AS repraesentiert und der Eintrag der  
	// jeweiligen zugehoerigen Gruppe (reduziertes Alphabet)
	// in diesem Alphabet
	// -> alphabet.cpp
	StringSet<String<AminoAcid>> Alphabete = GET_ALPHABET(comVal.numb_alp,comVal.size_alp);
	
	StringSet<String<Dna5>> Reads;
	StringSet<String<char>> ReadID;
	// fuellen der Container
	// -> own_functions.h 
	if (GET_DATAS(comVal.fastq_file,Reads,ReadID)) return 1;
	
	StrSetSA Proteine;
	StringSet<String<char>> ProteinID;
	// fuellen der Container
	// -> own_functions.h
	
	if (GET_DATAS(comVal.fasta_file,Proteine,ProteinID)) return 1;
	
	//
	if (FIND_MATCHES_FOR_ALL(Reads,ReadID,Proteine,ProteinID,comVal.seed,Alphabete)==1) return 1;
	*/
	return 0;
}
// MAIN FUNKTION -------------------------------------------------------------------------------