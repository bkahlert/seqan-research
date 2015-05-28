// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 19.APRIL.2013

#include "own_functions.h"

int hash(int x,int y,int z){
	int erg = (x*16)+(y*4)+z;
	return erg;
}

int get_AminoAcidPos(int pos){
	int x [64] = {4,15,4,15,17,17,17,17,5,14,5,14,11,11,0,11,16,6,16,6,8,8,8,8,5,5,5,5,10,10,10,10,2,3,2,3,7,7,7,7,12,12,12,12,9,9,9,9,20,18,20,18,14,14,14,14,20,1,19,1,10,13,10,13};
	return x[pos];
}

int translate(String<Dna> & aktual_triplet,String<int> & Alphabete){
	int hash_value = hash((int) aktual_triplet[0],(int) aktual_triplet[1],(int) aktual_triplet[2]);
	int gruppe = Alphabete[get_AminoAcidPos(hash_value)];
	cout<<aktual_triplet <<"\t"<<gruppe<<endl;

	return gruppe;
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
