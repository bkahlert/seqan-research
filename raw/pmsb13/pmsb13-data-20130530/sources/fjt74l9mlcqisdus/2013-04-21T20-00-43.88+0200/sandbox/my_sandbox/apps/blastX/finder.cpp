// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 19.APRIL.2013

#include "own_functions.h"

// 
void findMatches(StringSet<String<Dna>> & Reads,StringSet<String<char>> & ReadID,StringSet<String<AminoAcid>> & Proteine,
StringSet<String<char>> & ProteinID, int seed,StringSet<String<int>> & Alphabete){
	
	for(int alp=0;alp<length(Alphabete);++alp){
		StringSet<String<int>> trans_reads;
		for(int reading_frame=0;reading_frame<6;++reading_frame){
			translate_reads(trans_reads,Reads,Alphabete[alp],reading_frame);
		}
		for (int p=0;p<length(trans_reads);p++) cout << trans_reads[p]<<endl;

	}




}
