// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 19.APRIL.2013

#include "own_functions.h"

int hash(int x,int y,int z){
	int erg = (x*16)+(y*4)+z;
	return erg;
}

int get_AminoAcidPos(int pos){
	vector<int> x (1,2,3,4,5,19384);
	cout <<x[0]<<endl;
	return 1;
}

int translate(String<Dna> & aktual_triplet,String<int> & Alphabete){
	int hash_value = hash((int) aktual_triplet[0],(int) aktual_triplet[1],(int) aktual_triplet[2]);
	get_AminoAcidPos(1);
	return 0;
}



void translate_reads(StringSet<String<int>> & trans_reads,StringSet<String<Dna>> & Reads,String<int> & Alphabete,int frame){
	int numb_reads = length(Reads);
	StringSet<String<Dna>>revComplement;
	if ((frame/3)==1){
		revComplement=Reads;
		reverseComplement(revComplement);
		frame -= 3;
	}
	for (int i=0;i<numb_reads;++i){
		String<int>integer_code;
		for (int triplet=0+frame;triplet<length(Reads[i])&& triplet+3<length(Reads[i]);triplet+=3){
			String<Dna> aktual_triplet;
			if ((frame/3)==0) aktual_triplet = infix(Reads[i],triplet,triplet+3);
			else aktual_triplet = infix(revComplement[i],triplet,triplet+3);
			appendValue(integer_code,translate(aktual_triplet,Alphabete));
		}
		appendValue(trans_reads,integer_code);

	}


}
