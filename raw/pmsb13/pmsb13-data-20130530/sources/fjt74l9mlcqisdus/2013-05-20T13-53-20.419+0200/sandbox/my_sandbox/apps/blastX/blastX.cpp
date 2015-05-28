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
	if (CHECK_VALUES(comVal,reads)==1) return 1;

	//Erstellung eines reduzierten Alphabets
	// -> alphabet.cpp
	//StringSet<String<AminoAcid> > alphabets = GET_ALPHABET_FORCE();
	StringSet<String<AminoAcid> > alphabets = GET_ALPHABET(comVal.size_alp, comVal.numb_alp);
	// Alle Reads werden uebersetzt und in trans_reads gespeichert	
	// Alle Proteine der Datenbank werden uebersetzt und in trans_proteine gespeichert
	// -> translate.cpp
	for(unsigned alp=0;alp<length(alphabets);++alp){
		StrSetSA trans_reads;
		StrSetSA trans_proteine;

		if (TRANSLATE_READS(reads,trans_reads,alphabets[alp],1)) return 1;
		if (TRANSLATE_DATABASE(trans_proteine,proteine,alphabets[alp])) return 1;
	
		// Matches der Laenge seed werden gesucht und in seed_found gespeichert
		// -> finder.cpp
		Match_found seed_found;
		if (FIND_MATCHES(trans_proteine,trans_reads,comVal.seed,comVal.hamming_distance,seed_found)) return 1;
/*
//-------------------------------------------------------------------------------------------------------	

	//Proteinsequenzen erstellen	
	StrSetSA proteine;
	appendValue(proteine,"APLLQDP");
	appendValue(proteine,"LVIRVMVLIS");
	//Mit Bosum30-Matrix auf ein mit 10 Zeichen reduziertes Alphabet uebersetzt	
	StrSetSA trans_proteine;
	appendValue(trans_proteine, "NLRRDGL");
	appendValue(trans_proteine, "RNNDNQNRNC");
	StrSetSA trans_reads;
	appendValue(trans_reads, "CRNC");
	appendValue(trans_reads, "NLRN");	
	appendValue(trans_reads, "DNQQ");
	appendValue(trans_reads, "CKCC");	
	//Erlaubte missmatches 	
	unsigned distance = 1;	
	//Gespeicherte Informationen zu den matches
	Match_found seed_found;
	if (find_matches(trans_proteine,trans_reads,comVal.seed,comVal.hamming_distance,seed_found)) return 1;
	for(unsigned i=0;i<length(trans_proteine);++i){cout<<"trans_proteine:"<<trans_proteine[i]<<endl;}
	for(unsigned i=0;i<length(trans_reads);++i){cout<<"trans_reads:"<<trans_reads[i]<<endl;}
	cout<<"position_read:"<<seed_found.position_read<<endl;
	cout<<"begin_read:"<<seed_found.begin_read<<endl;
	cout<<"end_read:"<<seed_found.end_read<<endl;
	cout<<"position_protein:"<<seed_found.position_protein<<endl;
	cout<<"begin_protein:"<<seed_found.begin_protein<<endl;
	cout<<"end_protein:"<<seed_found.end_protein<<endl;
	Match_found verify_found;
	
	if(verify_all(seed_found, distance, trans_proteine, trans_reads, verify_found, proteine)) return 1;
	cout<<"position_read:"<<verify_found.position_read<<endl;
	cout<<"begin_read:"<<verify_found.begin_read<<endl;
	cout<<"end_read:"<<verify_found.end_read<<endl;
	cout<<"position_protein:"<<verify_found.position_protein<<endl;
	cout<<"begin_protein:"<<verify_found.begin_protein<<endl;
	cout<<"end_protein:"<<verify_found.end_protein<<endl;


//---------------------------------------------------------------------------------------------------*/		
		// Matches werden in eine Textdatei geschrieben
		// -> output.cpp
		WRITE_TO_FILE(seed_found, proteinID, readID,comVal.numb_alp,alp);
		
		// Matches werden verifiziert und in verify_found gespeichert
		// -> verify.cpp	
		Match_found verify_found;
		if (VERIFY_ALL(seed_found,comVal.hamming_distance,trans_proteine,trans_reads,verify_found,proteine)) return 1;
		
		// Verifizierte matches werden in eine Textdatei geschrieben
		// -> output.cpp
		//WRITE_TO_FILE(verify_found, proteinID, readID, reads,comVal.numb_alp,alp);
	}

	return 0;
}
// MAIN FUNKTION -------------------------------------------------------------------------------
