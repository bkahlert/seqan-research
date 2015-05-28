// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

#include "own_functions.h"

// 
void FIND_MATCHES_FOR_ALL(StringSet<String<Dna5>> & Reads,StringSet<String<char>> & ReadID,StrSetSA & Proteine,
StringSet<String<char>> & ProteinID, int seed,StrSetSA & Alphabete){
	// Schleife fuer jedes Alphabet 
	for(int alp=0;alp<length(Alphabete);++alp){
		StrSetSA trans_reads;
		// Schleife fuer jedes Read // Uebersetzung
		for (int read=0;read<length(Reads);++read){
			// Schleife fuer jedes moegliche Reading-Frame 
			for(int reading_frame=0;reading_frame<6;++reading_frame){
				// Funktion in translate.cpp
				translate_reads(trans_reads,Reads[read],Alphabete[alp],reading_frame);
			}
		}
		
		for (int p=0;p<length(trans_reads);p++) cout << trans_reads[p]<<endl;
		cout <<endl;

		StrSetSA trans_proteine;
		// Schleife fuer jedes Protein // Uebersetzung
		for (int protein=0;protein<length(Proteine);++protein){
			// Funktion in translate.cpp
			translate_datenbank(trans_proteine,Proteine[protein],Alphabete[alp]);
		}
		
		for (int p=0;p<length(trans_proteine);p++) cout << trans_proteine[p]<<endl;
		
		
		Index <StrSetSA> my_index(trans_proteine);
	
		Index <StrSetSA> my_index2(trans_reads);
		HammingDistance x = 1;
		
		//Finder <Index<StrSetSA>,Backtracking<HammingDistance>> finder(my_index);
		//Pattern <Index<StrSetSA>> my_pattern(my_index2);
		//while(find(finder,my_index2,0)){
		//	cout <<"["<<beginPosition(finder)<<";"<<endPosition(finder)<<")\t"<<infix(finder)<<endl;
		//}
		
	}
}




	
		