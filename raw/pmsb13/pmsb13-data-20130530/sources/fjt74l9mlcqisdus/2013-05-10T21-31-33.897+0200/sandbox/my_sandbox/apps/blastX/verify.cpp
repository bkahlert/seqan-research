// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 

#include "own_functions.h" 


Pair<unsigned int, unsigned int> get_position_in_prot(TAlign & align){
	Pair<unsigned int, unsigned int> return_pos;
	//for (int i=0;i<length())
	cout << row(alin,0)<<endl;
	cout << row(alin)<<endl;
	
	return return_pos;
	
}

int known_position(int & begin,int & end,StringSet<Pair<unsigned int, unsigned int>> & pair_position){
	for (int i=0;i<length(pair_position);++i){
		if (pair_position[i].i1 == begin && pair_position[i].i2 == end) return 0;
	}
	return 1;
}



int verify_seed_match(String<AminoAcid> protein,String<AminoAcid> & read,Match_found & verify_found){
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),protein);
	assignSource(row(align,1),read);
	int score = globalAlignment(align,Score<int,Simple>(1,-1,-1),NeedlemanWunsch());
	Pair<unsigned int, unsigned int> position = get_position_in_prot(align);
	
	
	if (score>0){
		cout << protein << "\t"<<read<<endl;	
		cout << "Score: " << score <<endl;
		cout << align <<endl;
	}
	return 0;

}	

	

int verify_all(Match_found & seed_found, int & distance,StrSetSA & trans_proteine,StrSetSA & trans_reads,Match_found & verify_found,StringSet<String<Dna> > & Reads,StrSetSA & Proteine){
	StrSetSA real_trans_reads;
	// alle Reads werden uebersetzt und in trans_reads gespeichert
	StringSet<String<AminoAcid> > forceAlphabet = GET_ALPHABET_FORCE();
	if (translate_reads(Reads,real_trans_reads,forceAlphabet[0],0)) return 1;
	int read_number = 0;
	do{
		StringSet<Pair<unsigned int, unsigned int>> pair_position;
		for (int i=0;i<length(seed_found.position_read);++i){
			if (read_number == seed_found.position_read[i]){
				String<AminoAcid> protein_areal;
				int begin = seed_found.begin_protein[i]-seed_found.begin_read[i];
				if (begin-distance>=0) begin-=distance;
				int end = seed_found.end_protein[i] + (length(real_trans_reads[read_number])-seed_found.end_read[i]);
				if (end+distance<=length(trans_proteine[seed_found.position_protein[i]])) end+=distance;
				if (begin>=0 && end<=length(trans_proteine[seed_found.position_protein[i]]) && known_position(begin,end,pair_position)){
					Pair<unsigned int,unsigned int> new_position (begin,end);
					appendValue(pair_position,new_position);
					verify_seed_match(infix(Proteine[seed_found.position_protein[i]],begin,end),real_trans_reads[read_number],verify_found);
				}
			}	
		}
		read_number ++;
	}while(read_number<length(trans_reads));
    
	return 0;
}
