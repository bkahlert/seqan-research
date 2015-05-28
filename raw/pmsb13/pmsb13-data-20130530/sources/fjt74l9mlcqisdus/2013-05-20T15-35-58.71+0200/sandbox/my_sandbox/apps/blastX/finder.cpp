// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

/**@brief finder.cpp
*Findet die Matches von den Reads mit der Proteindatenbank
*/

#include "own_functions.h"



int get_read_position(unsigned & pattern_pos,StringSet<unsigned> & found_reads,unsigned binary){
	/**@brief get_read_position Gibt entweder die Position des reads wieder oder die Startposition des seeds im reads
	*@param pattern_pos Position vom gematchten seed
	*@param found_reads Anzahl an seeds pro read
	*@param binary Bei 1 gibt die Funktion die Position des reads wieder und bei 0 die Position des seeds im read
	*@return Gibt entweder die Position des reads wieder oder die Startposition des seeds im reads
	*/
	unsigned sum = 0;
	for (unsigned count=0;count<length(found_reads);++count){
		unsigned sum_old = sum;
		sum+=found_reads[count];
		if (pattern_pos<sum && binary) return count;
		if (pattern_pos<sum && !binary){
			int read_begin = (sum-sum_old)-(sum-pattern_pos);
			return read_begin;
		}
		
	}
	return -1;
}

int append_to_match_found(Match_found & seed_found,Pair<unsigned int,unsigned int> & begin_found_prot,Pair<unsigned int,unsigned int> & end_found_prot,Pair<unsigned int,unsigned int> & pattern_found,StringSet<unsigned int> & found_reads){
	/**@brief append_to_match_found speichert Matchinformationen in der Klasse Match_found
	*@param seed_found Instanz der Klasse Match_found zur Speicherung der Matchinformationen der seeds
	*@param begin_found_prot Speichert die Anfangsposition des matches im Protein und die Position des Proteins selbst
	*@param end_found_prot Speichert die Endposition des matches im Protein und die Position des Proteins selbst
	*@param pattern_found Position des Patterns
	*@param found_reads Anzahl an seeds pro read
	*return Gibt boolschen Wert zurueck
	*/
	unsigned int protein_pos = begin_found_prot.i1;
	unsigned int protein_begin = begin_found_prot.i2;
	unsigned int protein_end = end_found_prot.i2;
	unsigned int pattern_pos = pattern_found.i1;
	int read_pos = get_read_position(pattern_pos,found_reads,1);
	int read_begin = get_read_position(pattern_pos,found_reads,0);

	if (read_pos==-1 || read_begin==-1) return 1;
	else{
	appendValue(seed_found.position_read,read_pos);
	appendValue(seed_found.begin_read,read_begin);
	appendValue(seed_found.end_read,read_begin + (protein_end-protein_begin));
	appendValue(seed_found.position_protein,protein_pos);
	appendValue(seed_found.begin_protein,protein_begin);
	appendValue(seed_found.end_protein,protein_end);
	}
	return 0;
}

int FIND_MATCHES(StrSetSA & trans_proteine,StrSetSA & trans_reads,unsigned & seed,unsigned & distance,Match_found & seed_found){
	/**@brief find_matches findet matches zwischen den Proteinen und seeds(verkuerzte Reads) 
	*@param trans_proteine Uebersetzte Proteine
	*@param trans_reads Uebersetzte Reads
	*@param seed Laenge der seeds
	*@param distance Anzahl an erlaubten missmatches
	*@param seed_found Speichert die gefundenen matches in der Klasse Match_found
	*@return Gibt boolschen Wert zurueck
	*/
	StrSetSA seed_reads;
	StringSet<unsigned> found_reads;
	for (unsigned read=0;read<length(trans_reads);++read){
		// nicht ganzes Pattern sondern Teilstücke werden gesucht
		unsigned begin_count = 0;
		for (unsigned begin=0;begin+seed<=length(trans_reads[read]);++begin){
			appendValue(seed_reads,infix(trans_reads[read],begin,begin+seed));
			begin_count++;
		}
		appendValue(found_reads,begin_count);
	}
	// Index ueber Proteindatenbank wird aufgebaut
	Index <StrSetSA, IndexSa<> > index_trans_proteine(trans_proteine);        
    Index <StrSetSA, IndexWotd<> > index_seed_reads(seed_reads);
    
    Finder <Index<StrSetSA, IndexSa<> >,Backtracking<HammingDistance> > finder_obj(index_trans_proteine);
    Pattern <Index<StrSetSA, IndexWotd<> >,Backtracking<HammingDistance> > pattern_obj(index_seed_reads,seed);   
	
	while(find(finder_obj,pattern_obj,distance)){
		//cout <<"end, begin:"<< beginPosition(finder_obj) <<"\t"<<endPosition(finder_obj)<<endl;
		Pair<unsigned int, unsigned int> begin_found_prot = beginPosition(finder_obj);
		Pair<unsigned int, unsigned int> end_found_prot = endPosition(finder_obj);
		Pair<unsigned int, unsigned int> pattern_found = position(pattern_obj);
		if (append_to_match_found(seed_found,begin_found_prot,end_found_prot,pattern_found,found_reads)) return 1;
	}
	return 0;
}

//int FIND_MATCHES_FOR_ALL(StringSet<String<Dna> > & reads,StringSet<String<char> > & readID,StrSetSA & proteine,StringSet<String<char> > & proteinID, int & seed,int & distance,StrSetSA & alphabets){
	/**@brief FIND_MATCHES_FOR_ALL uebersetzt, sucht matches, verifiziert und gibt die Ergebnisse in einer Textdatei aus 
	*@param reads Original Reads
	*@param readID IDs der original Reads
	*@param proteine Original Proteine aus der Proteindatenbank
	*@param proteinID IDs der Proteine 
	*@param seed Laenge der seeds
	*@param distance Anzahl an erlaubten missmatches
	*@param alphabet Reduziertes Alphabet
	*@return Gibt boolschen Wert zurueck
	*/
	// Schleife fuer jedes Alphabet 
/*	for(int alp=0;alp<length(alphabets);++alp){
		StrSetSA trans_reads;
		// alle Reads werden uebersetzt und in trans_reads gespeichert
		if (translate_reads(reads,trans_reads,alphabets[alp],1)) return 1;
	
		StrSetSA trans_proteine;
		// alle Proteine der Datenbank werden uebersetzt und in trans_proteine gespeichert 	
		if (translate_database(trans_proteine,proteine,alphabets[alp])) return 1;
	
		// Klasse angelegt und matches der laenge seed gesucht
		Match_found seed_found;
		if (find_matches(trans_proteine,trans_reads,seed,distance,seed_found)) return 1;
		
		write_to_file(seed_found, proteinID, readID);
		
		Match_found verify_found;
		if (verify_all(seed_found,distance,trans_proteine,trans_reads,verify_found,reads,proteine,alphabets)) return 1;
		
		write_to_file(verify_found, proteinID, readID, reads);
	}
	return 0;
}
*/




	
		
