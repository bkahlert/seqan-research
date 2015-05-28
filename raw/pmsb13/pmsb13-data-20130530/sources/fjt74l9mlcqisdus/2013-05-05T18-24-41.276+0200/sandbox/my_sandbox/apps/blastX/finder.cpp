// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

#include "own_functions.h"


// Speichert Match Informationen in Match_found klasse
void append_to_match_found(Match_found & seed_found,int read,int begin,int seed,Pair<unsigned int,unsigned int> & begin_found_prot,Pair<unsigned int,unsigned int> & end_found_prot){
	appendValue(seed_found.position_read,read);
	appendValue(seed_found.begin_read,begin);
	appendValue(seed_found.end_read,begin+seed);
	appendValue(seed_found.position_protein,begin_found_prot.i1);
	appendValue(seed_found.begin_protein,begin_found_prot.i2);
	appendValue(seed_found.end_protein,end_found_prot.i2);
}

// findet matches der laenge seed fuer ein datensatz von trans reads und einen datensatz trans proteine und speichert
// diese in eine Klasse von Match_found
int find_matches(StrSetSA & trans_proteine,StrSetSA & trans_reads,int seed,Match_found & seed_found){
	// Index ueber Proteindatenbank wird aufgebaut
	Index <StrSetSA> index_trans_proteine(trans_proteine);
	Finder <Index<StrSetSA>> finder_in_proteine(index_trans_proteine);
	// Fuer jedes Reads
	for (int read=0;read<length(trans_reads);++read){
		// nicht ganzes Pattern sondern Teilstücke werden gesucht
		for (int begin=0;begin+seed<=length(trans_reads[read]);++begin){
			while(find(finder_in_proteine,infix(trans_reads[read],begin,begin+seed))){
				Pair<unsigned int, unsigned int> begin_found_prot = beginPosition(finder_in_proteine);
				Pair<unsigned int, unsigned int> end_found_prot = endPosition(finder_in_proteine);
				append_to_match_found(seed_found,read,begin,seed,begin_found_prot,end_found_prot);
				
			}
			clear(finder_in_proteine);
		}
	}
	//for (int i=0;i<length(seed_found.position_read);++i) 
	//cout << seed_found.position_read << endl <<seed_found.begin_read << endl <<seed_found.end_read << endl <<
	//	seed_found.position_protein << endl <<seed_found.begin_protein << endl <<seed_found.end_protein << endl;
	return 0;
}

// fuer jedes alphabet wird übersetzt, seed matches gesucht, verifiziert und in einer textdatei ausgegeben
int FIND_MATCHES_FOR_ALL(StringSet<String<Dna>> & Reads,StringSet<String<char> > & ReadID,StrSetSA & Proteine,
StringSet<String<char> > & ProteinID, int seed,StrSetSA & Alphabete){
	// Schleife fuer jedes Alphabet 
	for(int alp=0;alp<length(Alphabete);++alp){
		StrSetSA trans_reads;
		// alle Reads werden uebersetzt und in trans_reads gespeichert
		if (translate_reads(Reads,trans_reads,Alphabete[alp])==1) return 1;
	
		StrSetSA trans_proteine;
		// alle Proteine der Datenbank werden uebersetzt und in trans_proteine gespeichert 	
		if (translate_database(trans_proteine,Proteine,Alphabete[alp])) return 1;
	
		// Klasse angelegt und matches der laenge seed gesucht
		Match_found seed_found;
		if (find_matches(trans_proteine,trans_reads,seed,seed_found)==1) return 1;
		
		void write_to_file(seed_found, ProteinID, ReadID);


	}
}




	
		
