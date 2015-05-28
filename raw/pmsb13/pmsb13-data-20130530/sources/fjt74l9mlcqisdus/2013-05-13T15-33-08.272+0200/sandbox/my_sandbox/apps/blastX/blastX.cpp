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
//------------------------------------------------------------------------------------------------
	for(unsigned i = 0; i < length(alphabets); ++i)
	{
		cout<<"Alphabet:"<<alphabets[i]<<endl;
	}

//------------------------------------------------------------------------------------------------
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
		write_to_file(seed_found, proteinID, readID);
		
		// Matches werden verifiziert und in verify_found gespeichert
		// -> verify.cpp	
		Match_found verify_found;
		if (verify_all(seed_found,comVal.hamming_distance,trans_proteine,trans_reads,verify_found,reads,proteine,alphabets)) return 1;
		
		// Verifizierte matches werden in eine Textdatei geschrieben
		// -> output.cpp
		write_to_file(verify_found, proteinID, readID, reads);

	}

	
	

//------------------------------------------------------------------------------------------------	
/*	for(unsigned i=0;i<length(trans_reads);++i){cout<<"trans_reads:"<<trans_reads[i]<<endl;}
	for(unsigned i=0;i<length(trans_proteine);++i){cout<<"trans_proteine:"<<trans_proteine[i]<<endl;}
	cout<<"position read:"<<seed_found.position_read<<endl;	
	cout<<"begin read:"<<seed_found.begin_read<<endl; 
	cout<<"end read:"<<seed_found.end_read<<endl;
	cout<<"positon protein:"<<seed_found.position_protein<<endl;
	cout<<"begin protein:"<<seed_found.begin_protein<<endl;
	cout<<"end protein:"<<seed_found.end_protein<<endl;
	cout<<"score:"<<seed_found.score<<endl;
*/	
//--------------------------------------------------------------------------------------------------
	// fuer alle alphabete werden matches gesucht verifiziert und ausgegeben
	/*if (FIND_MATCHES_FOR_ALL(Reads,ReadID,Proteine,ProteinID,comVal.seed,comVal.hamming_distance,alphabet)==1) return 1;*/
	
	return 0;
}
// MAIN FUNKTION -------------------------------------------------------------------------------
