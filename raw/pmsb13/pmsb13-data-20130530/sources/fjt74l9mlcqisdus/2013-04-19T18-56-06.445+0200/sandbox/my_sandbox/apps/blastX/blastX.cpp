// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 19.APRIL.2013

#include "own_functions.h"

// MAIN FUNKTION -------------------------------------------------------------------------------
int main(int argc, char const ** argv){
	// Klasse Variable enth�lt alle Parameter die der
	// Benutzer �ber Kommandozeile festlegen kann
	// -> own_functions.h
	Variable comVal;
	// comVal wird mit default Parameter initialisiert
	// -> parse_arguments.cpp
	dafault_values(comVal);
	// Eingabe vom Benutzer wird gelesen und default
	// Parameter wenn Eingabe vorhanden ist ueberschrieben
	if (argc>=2){
		// -> parse_arguments.cpp
		if (PARSE_ARGUMENTS(argc,argv,comVal))
			return 1;
	}
	// Alphabet enth�lt nach getAlphabet numb_alp mal 
	// String<int> der laenge zwanzig wobei jede Position
	// eine AS repraesentiert und der Eintrag der  
	// jeweiligen zugehoerigen Gruppe (reduziertes Alphabet)
	// in diesem Alphabet
	// -> alphabet.cpp
	StringSet<String<int>> Alphabete;
	getAlphabet(comVal.numb_alp,comVal.size_alp,Alphabete);
	
	StringSet<String<Dna>> Reads;
	StringSet<String<char>> ReadID;
	// fuellen der Container
	// -> own_functions.h 
	getDatas(comVal.fastq_file,Reads,ReadID);
	
	StringSet<String<AminoAcid>> Proteine;
	StringSet<String<char>> ProteinID;
	// fuellen der Container
	// -> own_functions.h
	getDatas(comVal.fasta_file,Proteine,ProteinID);
	
	//
	findMatches(Reads,ReadID,Proteine,ProteinID,comVal.seed,Alphabete);
	return 0;
}
// MAIN FUNKTION -------------------------------------------------------------------------------