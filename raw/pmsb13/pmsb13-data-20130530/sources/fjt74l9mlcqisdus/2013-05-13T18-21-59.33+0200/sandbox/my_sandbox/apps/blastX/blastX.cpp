// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

#include "own_functions.h"

// MAIN FUNKTION -------------------------------------------------------------------------------
int main(int argc, char const ** argv){
	
	// Klasse Variable enthält alle Parameter die der
	// Benutzer über Kommandozeile festlegen kann
	// -> own_functions.h
	Variable comVal;
	
	// comVal wird mit default Parameter initialisiert
	// -> parse_arguments.cpp
	DEFAULT_VALUES(comVal);
	
	// Eingabe vom Benutzer wird gelesen und default
	// Parameter wenn Eingabe vorhanden ist ueberschrieben
	// -> parse_arguments.cpp
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv,comVal))
			return 1;
	}
	
	// Fuellen der Read Container
	// -> own_functions.h
	StringSet<String<Dna> > reads;
	StringSet<String<char> > readID; 
	if (GET_DATAS(comVal.fastq_file,reads,readID)) return 1;

	// Fuellen der Protein Container
	// -> own_functions.h	
	StrSetSA proteine;
	StringSet<String<char> > proteinID;
	if (GET_DATAS(comVal.fasta_file,proteine,proteinID)) return 1;

	// Kontrolle der Eingabeparameter auf Machbarkeit
	// -> parse_arguments.cpp
	if (check_values(comVal,reads)==1) return 1;

	//Erstellung eines reduzierten Alphabets
	// -> alphabet.cpp
	//StringSet<String<AminoAcid> > alphabets = GET_ALPHABET_FORCE();
	StringSet<String<AminoAcid> > alphabets = GET_ALPHABET(comVal.size_alp, comVal.numb_alp);
	// Alle Reads werden uebersetzt und in trans_reads gespeichert	
	// Alle Proteine der Datenbank werden uebersetzt und in trans_proteine gespeichert
	// -> translate.cpp
	
	for(int alp=0;alp<length(alphabets);++alp){
		StrSetSA trans_reads;
		StrSetSA trans_proteine;

		if (translate_reads(reads,trans_reads,alphabets[alp],1)) return 1;
		if (translate_database(trans_proteine,proteine,alphabets[alp])) return 1;
	
		// Matches der Laenge seed werden gesucht und in seed_found gespeichert
		// -> finder.cpp
		Match_found seed_found;
		if (find_matches(trans_proteine,trans_reads,comVal.seed,comVal.hamming_distance,seed_found)) return 1;
		
		// Matches werden in eine Textdatei geschrieben
		// -> output.cpp
		write_to_file(seed_found, proteinID, readID,comVal.numb_alp,alp);
		
		// Matches werden verifiziert und in verify_found gespeichert
		// -> verify.cpp	
		Match_found verify_found;
		if (verify_all(seed_found,comVal.hamming_distance,trans_proteine,trans_reads,verify_found,reads,proteine,alphabets)) return 1;
		
		// Verifizierte matches werden in eine Textdatei geschrieben
		// -> output.cpp
		write_to_file(verify_found, proteinID, readID, reads,comVal.numb_alp,alp);

	}

	return 0;
}
// MAIN FUNKTION -------------------------------------------------------------------------------
