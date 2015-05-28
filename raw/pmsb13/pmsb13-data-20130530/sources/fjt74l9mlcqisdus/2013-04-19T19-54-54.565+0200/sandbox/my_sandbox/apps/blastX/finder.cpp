// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 19.APRIL.2013

#include "own_functions.h"

// 
void findMatches(StringSet<String<Dna>> & Reads,StringSet<String<char>> & ReadID,StringSet<String<AminoAcid>> & Proteine,
StringSet<String<char>> & ProteinID, int seed,StringSet<String<int>> & Alphabete){
	
	for(int alp=0;alp<length(Alphabete);++alp){
		for(int reading_frame=0;reading_frame<3;++reading_frame){
			StringSet<String<int>> trans_reads;
			translate_reads(trans_reads,Reads,Alphabete[alp],reading_frame);
			cout << length(trans_reads)<<endl;
			cout << trans_reads[0]<<endl;
			cout << trans_reads[1]<<endl;
			int g=0;
			cout <<g/3<<endl;
			g=4;
			cout<<g/3<<endl;
			g=1;
			cout <<g/3<<endl;
			g=5;
			cout<<g/3<<endl;
		}
	}




}
