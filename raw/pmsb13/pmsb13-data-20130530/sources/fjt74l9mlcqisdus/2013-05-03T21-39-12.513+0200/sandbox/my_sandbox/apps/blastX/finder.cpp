// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

#include "own_functions.h"
tringSet<String<char>> id_read; /**< ID Reads */
	StringSet<int> position_read; /**< Position Read in StringSet trans_reads */
	StringSet<int> begin_read; /**< begin in read */
	StringSet<int> end_read; /**< end in read */
	StringSet<String<char>> id_protein;
	StringSet<int> position_protein; /**< begin in protein */
	StringSet<int> begin_protein; /**< begin in protein */
	StringSet<int> end_protein; /**< end in protein */


//
int find_matches(StrSetSA & trans_proteine,StrSetSA & trans_reads,int seed,Match_found seed_found){
	Index <StrSetSA> index_trans_proteine(trans_proteine);
	Finder <Index<StrSetSA>> finder_in_proteine(index_trans_proteine);
	// Fuer jedes Reads
	for (int read=0;read<length(trans_reads);++read){
		for (int begin=0;begin+seed<=length(trans_reads[read]);begin+=seed){
			cout <<"infix "<<infix(trans_reads[read],begin,begin+seed)<<"\t"<<begin<<"\t"<<begin+seed<<endl;
			while(find(finder_in_proteine,infix(trans_reads[read],begin,begin+seed))){
				Pair<unsigned int, unsigned int> begin_found_prot = beginPosition(finder_in_proteine);
				Pair<unsigned int, unsigned int> end_found_prot = beginPosition(finder_in_proteine);
				//cout << x.i1<<endl;
				//Value<unsigned int,0>::Type y (x)
				appendValue(seed_found.position_read,read);
				appendValue(seed_found.begin_read,begin);
				appendValue(seed_found.end_read,begin);

				appendValue(seed_found.position_protein,begin_found_prot.i1);
				appendValue(seed_found.begin_protein,begin_found_prot.i2);
				appendValue(seed_found.end_read,end_found_prot.i2);
				//cout <<"["<<beginPosition(finder_in_proteine)<<";"<<endPosition(finder_in_proteine)<<")\t"<<
				//infix(finder_in_proteine)<<endl;
			}
			clear(finder_in_proteine);
		}
	}
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




	
		
