// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

/**@brief translate.cpp
*Beinhaltet Funktion, die sowohl die Reads als auch die Proteindatenbank uebersetzt.
*/

#include "own_functions.h"


int hash(int x,int y,int z){
	/**@brief hash gibt fuer jedes Codon einen eindeutigen integer wieder.
	*@param x Zahlenwert der Aminosäure an 1. Position 
	*@param y Zahlenwert der Aminosäure an 2. Position
	*@param z Zahlenwert der Aminosäure an 3. Position
	*@return erg gibt die Position in einem Array an , welcher die Aminosaeure gespeichert hat 
	*/
	int erg = (x*16)+(y*4)+z;
	if (erg>=0 && erg<=63) return erg;
	else {
		cerr << "invalid input in hash-function" <<endl;
		return -1;
	}
}

int get_Amino_Acid_Pos(int pos){
	/**@brief get_Amino_Acid_Pos beinhaltet ein Array, das die Adresse der jeweiligen Aminosaeure beinhaltet.
	*@param pos Position der gesuchten Aminosaeure
	*@return x[pos] Gibt die Adresse der Aminosaeure zurueck
	*/
	if (pos>=0 && pos<=63){
		int x [64] = {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,7,7,7,7,19,19,19,19,20,18,20,18,15,15,15,15,20,4,17,4,10,13,10,13};
		return x[pos];
	}
	else{
		cerr << "invalid input in get_Amino_Acid_Pos"<<endl;
		return -1;
	}
}


int get_translate_from_codon(String<Dna> & actual_triplet,String<AminoAcid> & alphabet,int reduced){
	/**@brief get_translate_from_codon gibt die Gruppennummer des reduzierten Alphabet zurueck der jeweiligen Amnosauere
	*@param actual_triplet Triplet, welches fuer die Aminosaeure codiert
	*@param alphabet Reduziertes Alphabet
	*@param reduced
	@return amino_pos Gruppennummer der Aminosaeure
	*/
	int hash_value = hash((int) actual_triplet[0],(int) actual_triplet[1],(int) actual_triplet[2]);
	if (hash_value!=-1){
		int amino_pos = get_Amino_Acid_Pos(hash_value);
		if (amino_pos!=-1 && reduced){
			int gruppe = alphabet[amino_pos];
			return gruppe;
		}
		if (amino_pos!=-1 && !reduced){
			return amino_pos;
		}
		else return -1;
	}
	else return -1;
}

int translate(StrSetSA & trans_reads,String<Dna> & read,String<AminoAcid> & alphabet,int frame,int reduced){
	/**@brief translate uebersetzt Read ueber die hash-Funktion in die jeweilige Aminosaeure
	*@param trans_reads Speichert die uebersetzten Reads
	*@param read Original Read
	*@param alphabet Reduziertes Alphabet
	*@param frame Moegliche reading Frames
	*@paramn reduced	
	*@retrun Gibt boolschen Wert zurueck
	*/
	String<Dna>revComplement;
	// frame = 0,1 oder 2 steht fuer forward 
	// frame = 3,4,oder 5 steht fuer reverse
	// complement gebildet und alle frames dafuer durchgegangen 
	int reverse = 1;
	if ((frame/3)==1){
		revComplement=read;
		reverseComplement(revComplement);
		frame -= 3;
		reverse = 0;
	}
	String<AminoAcid>integer_code;
	// Schleife durch Reads mit Frame Verschiebung
	for (int triplet=0+frame;triplet<length(read) && triplet+3<=length(read);triplet+=3){
		String<Dna> actual_triplet;
		if (reverse) actual_triplet = infix(read,triplet,triplet+3);
		else actual_triplet = infix(revComplement,triplet,triplet+3);
		int translate = get_translate_from_codon(actual_triplet,alphabet,reduced);
		if (translate!=-1) appendValue(integer_code,translate);
		else return 1;
	}
	appendValue(trans_reads,integer_code);
	return 0;
}

int translate_reads(StringSet<String<Dna> > & reads,StrSetSA & trans_reads,String<AminoAcid> & alphabet,int reduced){	
	/**@brief translate_reads uebersetzt Read fuer alle 6 reading Frames 
	*@param reads Original Reads
	*@param trans_reads Speichert die uebersetzten Reads
	*@param alphabet Reduziertes Alphabet
	*@param reduced
	@return Gibt boolschen Wert zurueck
	*/
	// Schleife fuer jedes Read // Uebersetzung
	for (int read=0;read<length(reads);++read){
		// Schleife fuer jedes moegliche Reading-Frame 
		for(int reading_frame=0;reading_frame<6;++reading_frame){
			// Funktion in translate.cpp
			if (translate(trans_reads,reads[read],alphabet,reading_frame,reduced)==1){
				cerr << "Programm fails"<<endl;
				return 1;
			}
		}
	}
	return 0;
}

int translate_database(StrSetSA & trans_proteine,StrSetSA & proteine,String<AminoAcid> & alphabet){	
	/**@brief translate_database uebersetzt die Proteine  in das vereinfachte Alphabet.
	*@param trans_proteine Speichert die uebersetzten Proteine 
	*@param proteine Proteine von der Datenbank
	*@param alphabet Reduziertes Alphabet
	*@return Gibt einen boolschen Wert zurueck
	*/
	// Schleife fuer jedes Protein // Uebersetzung
	for (int protein=0;protein<length(proteine);++protein){
		String<AminoAcid>new_amino_code;
		String<AminoAcid>old_amino_code= proteine[protein];
		for (int position=0;position<length(old_amino_code);++position){
			int amino = (int)old_amino_code[position];
			if (amino>=0 && amino<=20) appendValue(new_amino_code,alphabet[amino]);
			else{
				cerr << "invalid entry in database"<<endl;
				return 1;
			}
		}
		appendValue(trans_proteine,new_amino_code);
	}
	return 0;
}
