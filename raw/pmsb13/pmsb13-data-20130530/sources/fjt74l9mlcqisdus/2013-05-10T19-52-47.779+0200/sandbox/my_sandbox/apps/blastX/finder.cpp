// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

#include "own_functions.h"

//
int get_read_position(unsigned int & pattern_pos,StringSet<int> & found_reads,int binary){
	int sum = 0;
	for (int count=0;count<length(found_reads);++count){
		int sum_old = sum;
		sum+=found_reads[count];
		if (pattern_pos<sum && binary) return count;
		if (pattern_pos<sum && !binary){
			int read_begin = (sum-sum_old)-(sum-pattern_pos);
			return read_begin;
		}
		
	}
	return -1;
}


// Speichert Match Informationen in Match_found klasse
int append_to_match_found(Match_found & seed_found,Pair<unsigned int,unsigned int> & begin_found_prot,Pair<unsigned int,unsigned int> & end_found_prot,Pair<unsigned int,unsigned int> & pattern_found,int & seed,StringSet<int> & found_reads){
	unsigned int protein_pos = begin_found_prot.i1;
	unsigned int protein_begin = begin_found_prot.i2;
	unsigned int protein_end = end_found_prot.i2;
	unsigned int pattern_pos = pattern_found.i1;
	int read_pos = get_read_position(pattern_pos,found_reads,1);
	int read_begin = get_read_position(pattern_pos,found_reads,0);

	if (read_pos==-1 || read_begin==-1) return 1;
	else{
	appendValue(seed_found.position_read,read_pos);
	appendValue(seed_found.begin_read,read_begin);
	appendValue(seed_found.end_read,read_begin + (protein_end-protein_begin));
	appendValue(seed_found.position_protein,protein_pos);
	appendValue(seed_found.begin_protein,protein_begin);
	appendValue(seed_found.end_protein,protein_end);
	}
	return 0;
}

// findet matches der laenge seed fuer ein datensatz von trans reads und einen datensatz trans proteine und speichert
// diese in eine Klasse von Match_found
int find_matches(StrSetSA & trans_proteine,StrSetSA & trans_reads,int & seed,int & distance,Match_found & seed_found){
	StrSetSA seed_reads;
	StringSet<int> found_reads;
	for (int read=0;read<length(trans_reads);++read){
		// nicht ganzes Pattern sondern Teilstücke werden gesucht
		int begin = 0;
		for (begin;begin+seed<=length(trans_reads[read]);++begin){
			appendValue(seed_reads,infix(trans_reads[read],begin,begin+seed));
		}
		appendValue(found_reads,begin+1);
	}
	// Index ueber Proteindatenbank wird aufgebaut
	Index <StrSetSA, IndexSa<> > index_trans_proteine(trans_proteine);        
    Index <StrSetSA, IndexWotd<> > index_seed_reads(seed_reads);
    
    Finder <Index<StrSetSA, IndexSa<> >,Backtracking<HammingDistance> > finder_obj(index_trans_proteine);
    Pattern <Index<StrSetSA, IndexWotd<> >,Backtracking<HammingDistance> > pattern_obj(index_seed_reads,seed);   
	
	while(find(finder_obj,pattern_obj,distance)){
		cout << beginPosition(finder_obj) <<"\t"<<endPosition(finder_obj)<<endl;
		Pair<unsigned int, unsigned int> begin_found_prot = beginPosition(finder_obj);
		Pair<unsigned int, unsigned int> end_found_prot = endPosition(finder_obj);
		Pair<unsigned int, unsigned int> pattern_found = position(pattern_obj);
		if (append_to_match_found(seed_found,begin_found_prot,end_found_prot,pattern_found,seed,found_reads)) return 1;
	}
	return 0;
}

// fuer jedes alphabet wird übersetzt, seed matches gesucht, verifiziert und in einer textdatei ausgegeben
int FIND_MATCHES_FOR_ALL(StringSet<String<Dna> > & Reads,StringSet<String<char> > & ReadID,StrSetSA & Proteine,
StringSet<String<char> > & ProteinID, int & seed,int & distance,StrSetSA & Alphabete){
	// Schleife fuer jedes Alphabet 
	for(int alp=0;alp<length(Alphabete);++alp){
		StrSetSA trans_reads;
		// alle Reads werden uebersetzt und in trans_reads gespeichert
		if (translate_reads(Reads,trans_reads,Alphabete[alp],1)) return 1;
	
		StrSetSA trans_proteine;
		// alle Proteine der Datenbank werden uebersetzt und in trans_proteine gespeichert 	
		if (translate_database(trans_proteine,Proteine,Alphabete[alp])) return 1;
	
		// Klasse angelegt und matches der laenge seed gesucht
		Match_found seed_found;
		if (find_matches(trans_proteine,trans_reads,seed,distance,seed_found)) return 1;
		
		write_to_file(seed_found, ProteinID, ReadID);
		
		Match_found verify_found;
		if (verify_all(seed_found,distance,trans_proteine,trans_reads,verify_found,Reads,Proteine)) return 1;

	}
	return 0;
}




	
		
