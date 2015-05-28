// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

#include "own_functions.h"

// 
void findMatches(StringSet<String<Dna>> & Reads,StringSet<String<char>> & ReadID,StringSet<String<AminoAcid>> & Proteine,
StringSet<String<char>> & ProteinID, int seed,StringSet<String<int>> & Alphabete){
	// Schleife fuer jedes Alphabet 
	for(int alp=0;alp<length(Alphabete);++alp){
		StringSet<String<int>> trans_reads;
		// Schleife fuer jedes Read // Uebersetzung
		for (int read=0;read<length(Reads);++read){
			// Schleife fuer jedes moegliche Reading-Frame 
			for(int reading_frame=0;reading_frame<6;++reading_frame){
				// Funktion in translate.cpp
				translate_reads(trans_reads,Reads[read],Alphabete[alp],reading_frame);
			}
		}
		for (int p=0;p<length(trans_reads);p++) cout << trans_reads[p]<<endl;
		// Schleife fuer jedes Protein // Uebersetzung
		StringSet<String<int>> trans_databank;
		String<AminoAcid> x = Proteine[0];
		for (int p=0;p<length(x);++p) cout << (int)x[p]<<"\t";
		//for (int protein=0;protein<length(Proteine);++protein){
			
		//}
	}
}
