// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 

#include "own_functions.h" 


int append_to_match_found(Match_found & verify_found, Pair<int,int> & intervall_protein,int 
position_protein,int & position_read,int score){
	/**@brief appen_to_match_found speichert match Informationen in die Match_found Klasse.
	*@param verify_found Instanz der Klasse Match_found speichert die Daten
	*@param intervall_protein
	*@param position_protein Positionen der gefundenen matches im Protein
	*@param position_read Postionen der gefundenen matches im Read
	*@param score Score des verifizierten matches	
	*@return Beim erfolgreichen Speichern der Daten in die Klasse wird eine 0 zurueckgegeben
	*/
	appendValue(verify_found.position_protein,position_protein);
	appendValue(verify_found.begin_protein,intervall_protein.i1);
	appendValue(verify_found.end_protein,intervall_protein.i2);
	appendValue(verify_found.position_read,position_read);
	appendValue(verify_found.score,score);
	return 0;
}

Pair<int,int> get_position_in_prot(TAlign & align,int & begin,int & end){
	/**@brief get_position_in_prot
	*@param align Ein Alignment
	*@param begin Anfangsposition des Alignments
	*@param end Endposition des Alignments
	*@return return_pos 
	*/
	Pair<int,int> return_pos;
	for (int i=0;i<length(row(align,1));++i){
		if (row(align,1)[i]!='-') {
			return_pos.i1 = begin+i;
			break;
		}
		if (i==length(row(align,1))-1) return_pos.i1 = -1;
	}
	for (int i=0;i<length(row(align,1));++i){
		if (row(align,1)[length(row(align,1))-i-1]!='-') {
			return_pos.i2 = end-i;
			break;
		}
		if (i==length(row(align,1))-1) return_pos.i2 = -1;

	}		
	return return_pos;
}

int known_position(int & begin,int & end,StringSet<Pair<unsigned int, unsigned int> > & pair_position){
	/**@brief known_position Ueberprueft ob die Anfangs- und Endpositionen bereits betrachtet wurden
	*@param begin Aktuelle Anfangsposition im Protein
	*@param end Aktuelle Endposition im Protein
	*@param pair_position Speichert bereits betrachtete Anfangs- und Endpositionen im Protein 
	*@return return_pos Gibt 1 zuruek wenn die Postionen noch nicht betrachtet wurden
	*/
	for (int i=0;i<length(pair_position);++i){
		if (pair_position[i].i1 == begin && pair_position[i].i2 == end) return 0;
	}
	return 1;
}

int verify_seed_match(String<AminoAcid> protein,String<AminoAcid> & read,Match_found & verify_found,int & begin, int & end, int position_protein,int & position_read,int distance){
	/**@brief verify_seed_match
	*@param protein
	*@param read
	*@param verify_found
	*@param begin
	*@param end
	*@param position_protein
	*@param position_read
	*@return Gibt einen boolschen Wert zurueck
	*/
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),protein);
	assignSource(row(align,1),read);
	int score = globalAlignment(align,Score<int,Simple>(1,-1,-1),NeedlemanWunsch());
	int untere_grenze = length(read)-(3*distance);
	if (score>untere_grenze){
		cout << protein << "\t"<<read<<endl;	
		cout << "Score: " << score <<endl;
		cout << align <<endl;
		Pair<int,int> intervall = get_position_in_prot(align,begin,end);
		if (intervall.i1!=-1 && intervall.i2!=-1) append_to_match_found(verify_found,intervall,position_protein,position_read,score);
		else return 1;
	}
	return 0;
}	

int verify_all(Match_found & seed_found, int & distance,StrSetSA & trans_proteine,StrSetSA & trans_reads,Match_found & verify_found,StrSetSA & alphabets){
	/**@biref verify_all Verifiziert die gefundenen matches und speichert diese
	*@param seed_found Enthaelt die bereits gefundenen matches
	*@param distance Anzahl an erlaubten Fehlern
	*@param trans_proteine Uebersetzte Proteinsequenzen
	*@param trans_reads Uebersetzte Reads
	*@param verify_found Speichert die verifizierten matches
	*@param reads Originalreads
	*@param proteine Originalproteinsequenzen
	*@param alphabets Reduziertes Alphabet
	*@return Beim erfolgreichen Verifizieren gibt es eine 0 zurueck
	*/
	int read_number = 0;
	do{
		StringSet<Pair<unsigned int, unsigned int> > pair_position;
		for (int i=0;i<length(seed_found.position_read);++i){
			if (read_number == seed_found.position_read[i]){
				//String<AminoAcid> protein_areal;
				int begin = seed_found.begin_protein[i]-seed_found.begin_read[i];
				if (begin-distance>=0) begin-=distance;
				int end = seed_found.end_protein[i] + (length(trans_reads[read_number])-seed_found.end_read[i]);
				if (end+distance<=length(trans_proteine[seed_found.position_protein[i]])) end+=distance;
				if (begin>=0 && end<=length(trans_proteine[seed_found.position_protein[i]]) && known_position(begin,end,pair_position)){
					Pair<unsigned int,unsigned int> new_position (begin,end);
					appendValue(pair_position,new_position);
					verify_seed_match(infix(trans_proteine[seed_found.position_protein[i]],begin,end),trans_reads[read_number],verify_found,begin,end,seed_found.position_protein[i],read_number,distance);
				}
			}	
		}
		read_number ++;
	}while(read_number<length(trans_reads));
    
	return 0;
}
