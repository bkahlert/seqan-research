// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

#include "own_functions.h"

// HASH-FUNKTION GIBT FUER JEDES CODON EIN EINDEUTIGEN INT WIEDER
// DER INT GIBT DIE POSITION IN EINEM ARRAY AN, WELCHER INDIREKT DIE
// AMINOSAEURE GESPEICHERT HAELT

int hash(int x,int y,int z){
	/**@brief Hashfunktion 
	*@param x int zahlenwert der Aminosäure an 1. Position 
	*@param y int zahlenwert der Aminosäure an 2. Position
	*@param z int zahlenwert der Aminosäure an 3. Position
	*@return int [0,64] 
	*/
	int erg = (x*16)+(y*4)+z;
	if (erg>=0 && erg<=63) return erg;
	else {
		cerr << "invalid input in hash-function" <<endl;
		return -1;
	}
}

// FUNKTION BEINHALTET EIN ARRAY MIT 64 POSITIONEN WELCHE DIE ADRESSE
// EINER AMINOSAEURE BEINHALTEN 
int get_Amino_Acid_Pos(int pos){
	if (pos>=0 && pos<=63){
		int x [64] = {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,7,7,7,7,19,19,19,19,20,18,20,18,15,15,15,15,20,4,17,4,10,13,10,13};
		return x[pos];
	}
	else{
		cerr << "invalid input in get-amino-acid-pos"<<endl;
		return -1;
	}
}

// BEKOMMT EIN CODON UND GIBT GRUPPENNUMMER DER JEWEILIGEN AMINOSÄURE ZURUECK
int get_translate_from_codon(String<Dna> & aktual_triplet,String<AminoAcid> & Alphabete,int reduziert){
	int hash_value = hash((int) aktual_triplet[0],(int) aktual_triplet[1],(int) aktual_triplet[2]);
	if (hash_value!=-1){
		int Amino_pos = get_Amino_Acid_Pos(hash_value);
		if (Amino_pos!=-1 && reduziert){
			int gruppe = Alphabete[Amino_pos];
			return gruppe;
		}
		if (Amino_pos!=-1 && !reduziert){
			return Amino_pos;
		}
		else return -1;
	}
	else return -1;
}

// BEKOMMT EIN READ UND UEBERSETZT ES ZUERST UEBER EINE HASH-FUNKTION IN DIE JEWEILIGE AMINOSAEURE
// DA READING FRAME NICHT BEKANNT IST AUCH IN ALLE 6 MOEGLICHEN READING FRAMES (WIRD NICHT GESPEICHERT)
// SONDERN NUR MIT DER GRUPPEN NUMMER DES JEWEILIGEN ALPHABETES IN TRANS_READS GESPEICHERT 
int translate(StrSetSA & trans_reads,String<Dna> & read,String<AminoAcid> & Alphabete,int frame,int reduziert){
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
		String<Dna> aktual_triplet;
		if (reverse) aktual_triplet = infix(read,triplet,triplet+3);
		else aktual_triplet = infix(revComplement,triplet,triplet+3);
		int translate = get_translate_from_codon(aktual_triplet,Alphabete,reduziert);
		if (translate!=-1) appendValue(integer_code,translate);
		else return 1;
	}
	appendValue(trans_reads,integer_code);
	return 0;
}

// Uebersetz die Reads 1:6
int translate_reads(StringSet<String<Dna> > & Reads,StrSetSA & trans_reads,String<AminoAcid> & Alphabete,int reduziert){	
// Schleife fuer jedes Read // Uebersetzung
	for (int read=0;read<length(Reads);++read){
		// Schleife fuer jedes moegliche Reading-Frame 
		for(int reading_frame=0;reading_frame<6;++reading_frame){
			// Funktion in translate.cpp
			if (translate(trans_reads,Reads[read],Alphabete,reading_frame,reduziert)==1){
				cerr << "Programm fails"<<endl;
				return 1;
			}
		}
	}
	return 0;
}

// BEKOMMT ALLE PROTEINSEQUENZEN UND UEBERSETZT DIESE IN DAS VEREINFACHTE ALPHABET
int translate_database(StrSetSA & trans_proteine,StrSetSA & Proteine,String<AminoAcid> & Alphabete){
	// Schleife fuer jedes Protein // Uebersetzung
	for (int protein=0;protein<length(Proteine);++protein){
		String<AminoAcid>new_amino_code;
		String<AminoAcid>old_amino_code= Proteine[protein];
		for (int position=0;position<length(old_amino_code);++position){
			int Amino = (int)old_amino_code[position];
			if (Amino>=0 && Amino<=20) appendValue(new_amino_code,Alphabete[Amino]);
			else{
				cerr << "invalid entry in database"<<endl;
				return 1;
			}
		}
		appendValue(trans_proteine,new_amino_code);
	}
	return 0;
}
