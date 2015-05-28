// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

#include "own_functions.h"

// MAIN FUNKTION -------------------------------------------------------------------------------
// ALLES FUNKTIONEN DIE VON DER MAIN AUFGERUFEN WERDEN SIND GROSS GESCHRIEBEN
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
	cout<<"GET DATAS"<<endl;
	StringSet<String<Dna> > reads;
	StringSet<String<char> > readID; 
	if (GET_DATAS(comVal.fastq_file,reads,readID)) return 1;

	// Fuellen der Protein Container
	// -> own_functions.h	
	StrSetSA proteine;
	StringSet<String<char> > proteinID;
	if (GET_DATAS(comVal.fasta_file,proteine,proteinID)) return 1;
	int database_length = length(proteine);

	// Kontrolle der Eingabeparameter auf Machbarkeit
	// -> parse_arguments.cpp
	cout<<"CHECK VALUES"<<endl;
	if (CHECK_VALUES(comVal,reads)==1) return 1;

	//Erstellung eines reduzierten Alphabets
	// -> alphabet.cpp
	cout<<"CREATE ALPHABET"<<endl;
	//StringSet<String<AminoAcid> > alphabets = GET_ALPHABET_FORCE();
	StringSet<String<AminoAcid> > alphabets = GET_ALPHABET(comVal.size_alp, comVal.numb_alp);
	
	// Zeitmessung
	clock_t start_time, end_time;
	start_time = clock();
    
	for(unsigned alp=0;alp<length(alphabets);++alp){
		// Alle Reads werden uebersetzt und in trans_reads gespeichert	
		// Alle Proteine der Datenbank werden uebersetzt und in trans_proteine gespeichert
		// -> translate.cpp		
		StrSetSA trans_reads;
		StrSetSA trans_proteine;
		cout << "TRANSLATE_READS_"<<alp<<endl;		
		if (TRANSLATE_READS(reads,trans_reads,alphabets[alp],1)) return 1;
		cout << "TRANSLATE_DATABASE_"<<alp<<endl;
		if (TRANSLATE_DATABASE(trans_proteine,proteine,alphabets[alp])) return 1;

		
		// Matches der Laenge seed werden gesucht und in seed_found gespeichert
		// -> finder.cpp
		cout<<"FIND MATCHES "<<alp<<endl;
		Match_found seed_found;
		if (FIND_MATCHES(trans_proteine,trans_reads,comVal.seed,comVal.hamming_distance,seed_found)) return 1;

		// Matches werden in eine Textdatei geschrieben
		// -> output.cpp
		//WRITE_TO_FILE(seed_found, proteinID, readID,comVal.numb_alp,alp);
		
		// Matches werden verifiziert und in verify_found gespeichert
		// -> verify.cpp		
		
		StrSetSA real_trans_reads;
		if (TRANSLATE_READS(reads,real_trans_reads,alphabets[alp],0)) return 1;
		
		cout << "VERIFY_ALL_"<<alp<<endl;
		Match_found verify_found;
		if (VERIFY_ALL(seed_found,comVal.hamming_distance,trans_proteine,real_trans_reads,verify_found,proteine,database_length)) return 1;
		
		
		// Verifizierte matches werden in eine Textdatei geschrieben
		// -> output.cpp
		cout<<"WRITE TO FILE "<<alp<<endl;
		WRITE_TO_FILE(verify_found, proteinID, readID, reads,comVal.numb_alp,alp);
	}
	cout << "RUN_COMPLETE"<<endl;	
	end_time = clock();
	cout << start_time << endl;
	cout << end_time << ':' << CLOCKS_PER_SEC << endl;
	cout <<"READY"<<endl;
	return 0;
}
// MAIN FUNKTION -------------------------------------------------------------------------------
