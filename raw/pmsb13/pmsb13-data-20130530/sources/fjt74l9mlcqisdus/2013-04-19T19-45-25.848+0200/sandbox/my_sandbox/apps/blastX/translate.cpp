// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 19.APRIL.2013

#include "own_functions.h"

int translate(String<Dna> & aktual_triplet,String<int> & Alphabete){
	return 0;
}



void translate_reads(StringSet<String<int>> & trans_reads,StringSet<String<Dna>> & Reads,String<int> & Alphabete,int frame){
	int numb_reads = length(Reads);
	for (int i=0;i<numb_reads;++i){
		String<int>integer_code;
		for (int triplet=0+frame;triplet<length(Reads[i])&& triplet+3<length(Reads[i]);triplet+=3){
			String<Dna> aktual_triplet = infix(Reads[i],triplet,triplet+3);
			appendValue(integer_code,translate(aktual_triplet,Alphabete));
		}
		appendValue(trans_reads,integer_code);
	}


}
