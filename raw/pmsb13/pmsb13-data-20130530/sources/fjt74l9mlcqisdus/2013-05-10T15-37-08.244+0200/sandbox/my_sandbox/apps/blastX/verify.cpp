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
int verify_seed_match(String<AminoAcid> & protein,String<AminoAcid> & read,Match_found & verify_found){
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),protein);
	assignSource(row(align,1),read);
	int score = globalAlignment(align,Score<int,Simple>(1,-1,-1),AlignConfig<true, true, true, true>(),MyerBitVector());
	if (score>???){
		get_position_in_prot(align);
		append_to_match_found();
	}
	cout << "Score: " << score <<endl;
	cout << align <<endl;
	return 0;

}	

	

int verify_all(Match_found & seed_found, int & distance,StrSetSA & trans_proteine,StrSetSA & trans_reads,Match_found & verify_found,int & read_length){
	
	for (int match = 0;match<length(seed_found.position_read);++match){
		String<AminoAcid> protein_areal;
		int begin = seed_found.begin_protein[match]-seed_found.begin_read[match]+distance;	//distance <=0
		int end = seed_found.end_protein + (read_length-seed_found.end_read)-distance;
		//verify_seed_match(infix(trans_protein[seed_found.position_protein[match]],begin,end),trans_reads[seed_found.position_read[match]],verify_found);

	
	
	}
	    

    return 0;
}
