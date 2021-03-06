// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

#include "own_functions.h"

// 
int FIND_MATCHES_FOR_ALL(StringSet<String<Dna5> > & Reads,StringSet<String<char> > & ReadID,StrSetSA & Proteine,
StringSet<String<char> > & ProteinID, int seed,StrSetSA & Alphabete){
	// Schleife fuer jedes Alphabet 
	for(int alp=0;alp<length(Alphabete);++alp){
		StrSetSA trans_reads;
		// Schleife fuer jedes Read // Uebersetzung
		for (int read=0;read<length(Reads);++read){
			// Schleife fuer jedes moegliche Reading-Frame 
			for(int reading_frame=0;reading_frame<6;++reading_frame){
				// Funktion in translate.cpp
				if (translate_reads(trans_reads,Reads[read],Alphabete[alp],reading_frame)==1){
					cerr << "Programm fails"<<endl;
					return 1;
				}
			}
		}
	
		for (int p=0;p<length(trans_reads);p++) cout <<trans_reads[p]<<endl;
		cout <<endl;

		StrSetSA trans_proteine;
		// Schleife fuer jedes Protein // Uebersetzung
		for (int protein=0;protein<length(Proteine);++protein){
			// Funktion in translate.cpp
			if (translate_datenbank(trans_proteine,Proteine[protein],Alphabete[alp])) return 1;
		}
		
		for (int p=0;p<length(trans_proteine);p++) cout << trans_proteine[p]<<endl;
		
	
		Index <StrSetSA> index_trans_proteine(trans_proteine);
		Finder <Index<StrSetSA>> finder_in_proteine(index_trans_proteine);

		for (int read=0;read<length(trans_reads);++read){
			for (int begin=0;begin+seed<=length(trans_reads[read]);begin+=seed){
				cout <<"infix "<<infix(trans_reads[read],begin,begin+seed)<<endl;
				while(find(finder_in_proteine,infix(trans_reads[read],begin,begin+seed))){
					cout <<"["<<beginPosition(finder_in_proteine)<<";"<<endPosition(finder_in_proteine)<<")\t"<<
						infix(finder_in_proteine)<<endl;
				}
				clear(finder_in_proteine);
			}
			
		}
		
		/*
		
		Index <StrSetSA> index_trans_proteine(trans_proteine);		
		Index <StrSetSA> index_trans_reads(x);
	
		Finder <Index<StrSetSA>,Backtracking<HammingDistance> > finder_obj(index_trans_proteine);
		Pattern <Index<StrSetSA>,Backtracking<HammingDistance> > pattern_obj(index_trans_reads,length(x));
	    
		while(find(finder_obj,pattern_obj,0)){
			cout <<"["<<beginPosition(finder_obj)<<";"<<endPosition(finder_obj)<<")\t"<<infix(finder_obj)<<endl;
		}
		*/
				
	}
}




	
		
