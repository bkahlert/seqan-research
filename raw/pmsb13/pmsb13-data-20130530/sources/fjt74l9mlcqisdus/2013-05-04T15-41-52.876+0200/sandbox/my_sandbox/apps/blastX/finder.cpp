// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

#include "own_functions.h"


//
void append_to_match_found(int read,int begin,int seed,Pair<unsigned int,unsigned int> & begin_found_prot,Pair<unsigned int,unsigned int> & end_found_prot){
	appendValue(seed_found.position_read,read);
	appendValue(seed_found.begin_read,begin);
	appendValue(seed_found.end_read,begin+seed);
	appendValue(seed_found.position_protein,begin_found_prot.i1);
	appendValue(seed_found.begin_protein,begin_found_prot.i2);
	appendValue(seed_found.end_protein,end_found_prot.i2);
}
//
int find_matches(StrSetSA & trans_proteine,StrSetSA & trans_reads,int seed,Match_found seed_found){
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
				append_to_match_found(read,begin,seed,begin_found_prot,end_found_prot);
				
			}
			clear(finder_in_proteine);
		}
	}
	//for (int i=0;i<length(seed_found.position_read);++i) 
	//cout << seed_found.position_read << endl <<seed_found.begin_read << endl <<seed_found.end_read << endl <<
	//	seed_found.position_protein << endl <<seed_found.begin_protein << endl <<seed_found.end_protein << endl;
	return 0;
}

// 
int FIND_MATCHES_FOR_ALL(StringSet<String<Dna5>> & Reads,StringSet<String<char> > & ReadID,StrSetSA & Proteine,
StringSet<String<char> > & ProteinID, int seed,StrSetSA & Alphabete){
	// Schleife fuer jedes Alphabet 
	for(int alp=0;alp<length(Alphabete);++alp){
		StrSetSA trans_reads;
		// alle Reads werden uebersetzt und in trans_reads gespeichert
		if (get_translate_for_all(Reads,trans_reads,Alphabete[alp])==1) return 1;
		// Ausgabe der Reads	
		for (int p=0;p<length(trans_reads);p++) cout <<trans_reads[p]<<endl;
		cout <<endl;
	
		StrSetSA trans_proteine;
		// alle Proteine der Datenbank werden uebersetzt und in trans_proteine gespeichert 	
		if (translate_datenbank(trans_proteine,Proteine,Alphabete[alp])) return 1;
		// Ausgabe der Proteine	
		for (int p=0;p<length(trans_proteine);p++) cout << trans_proteine[p]<<endl;
		cout <<endl;
	
		Match_found seed_found;
		if (find_matches(trans_proteine,trans_reads,seed,seed_found)==1) return 1;
		
		
		/*
		
		Index <StrSetSA> index_trans_proteine(trans_proteine);		
		Index <StrSetSA> index_trans_reads(x);
	
		Finder <Index<StrSetSA>,Backtracking<HammingDistance> > finder_obj(index_trans_proteine);
		Pattern <Index<StrSetSA>,Backtracking<HammingDistance> > pattern_obj(index_trans_reads,length(x));
	    
		while(find(finder_obj,pattern_obj,0)){
			cout <<"["<<beginPosition(finder_obj)<<";"<<endPosition(finder_obj)<<")\t"<<infix(finder_obj)<<endl;
		}
		*/
				
	}
}




	
		
