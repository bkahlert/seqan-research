// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 

#include "own_functions.h" 


/*

	get_position_in_prot(align){
	for (int pos=0;pos<=length(row(align));++pos){
		row(alin,pos);
	}
	
}
*/	
	

int verify_seed_match(String<AminoAcid> protein,String<AminoAcid> & read,Match_found & verify_found){
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),protein);
	assignSource(row(align,1),read);
	int score = globalAlignment(align,Score<int,Simple>(1,-1,-1),NeedlemanWunsch());
	/*get_position_in_prot(align);
	append_to_match_found();
	*/
	if (score>1){
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
	
	for (int match = 0;match<length(seed_found.position_read);++match){
		String<AminoAcid> protein_areal;
		int begin = seed_found.begin_protein[match]-seed_found.begin_read[match]+distance;	//distance <=0
		int end = seed_found.end_protein[match] + (length(trans_reads[seed_found.position_read[match]])-seed_found.end_read[match])-distance;
		
		if (begin>=0 && end<=length(trans_proteine[seed_found.position_protein[match]])){
			verify_seed_match(infix(Proteine[seed_found.position_protein[match]],begin,end),real_trans_reads[seed_found.position_read[match]],verify_found);
		}
	
	
	}
	
    return 0;
}
