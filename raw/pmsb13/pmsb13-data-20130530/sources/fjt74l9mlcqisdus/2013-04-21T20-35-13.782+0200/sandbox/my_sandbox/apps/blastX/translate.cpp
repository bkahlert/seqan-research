// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 19.APRIL.2013

#include "own_functions.h"

// HASH-FUNKTION GIBT FUER JEDES CODON EIN EINDEUTIGEN INT WIEDER
// DER INT GIBT DIE POSITION IN EINEM ARRAY AN, WELCHER INDIREKT DIE
// AMINOSAEURE GESPEICHERT HAELT
int hash(int x,int y,int z){
	int erg = (x*16)+(y*4)+z;
	return erg;
}

// FUNKTION BEINHALTET EIN ARRAY MIT 64 POSITIONEN WELCHE DIE ADRESSE
// EINER AMINOSAEURE BEINHALTEN 
int get_AminoAcidPos(int pos){
	int x [64] = {4,15,4,15,17,17,17,17,5,14,5,14,11,11,0,11,16,6,16,6,8,8,8,8,5,5,5,5,10,10,10,10,2,3,2,3,7,7,7,7,12,12,12,12,9,9,9,9,20,18,20,18,14,14,14,14,20,1,19,1,10,13,10,13};
	return x[pos];
}

// BEKOMMT EIN CODON UND GIBT GRUPPENNUMMER DER JEWEILIGEN AMINOSÄURE ZURUECK
int translate(String<Dna> & aktual_triplet,String<int> & Alphabete){
	int hash_value = hash((int) aktual_triplet[0],(int) aktual_triplet[1],(int) aktual_triplet[2]);
	int gruppe = Alphabete[get_AminoAcidPos(hash_value)];
	return gruppe;
}

// BEKOMMT EIN READ UND UEBERSETZT ES ZUERST UEBER EINE HASH-FUNKTION IN DIE JEWEILIGE AMINOSAEURE
// DA READING FRAME NICHT BEKANNT IST AUCH IN ALLE 6 MOEGLICHEN READING FRAMES (WIRD NICHT GESPEICHERT)
// SONDERN NUR MIT DER GRUPPEN NUMMER DES JEWEILIGEN ALPHABETES IN TRANS_READS GESPEICHERT 
void translate_reads(StringSet<String<int>> & trans_reads,String<Dna> & read,String<int> & Alphabete,int frame){
	String<Dna>revComplement;
	if ((frame/3)==1){
		revComplement=read;
		reverseComplement(revComplement);
		frame -= 3;
	}
	String<int>integer_code;
	for (int triplet=0+frame;triplet<length(read) && triplet+3<length(read);triplet+=3){
		String<Dna> aktual_triplet;
		if ((frame/3)==0) aktual_triplet = infix(read,triplet,triplet+3);
		else aktual_triplet = infix(revComplement,triplet,triplet+3);
		appendValue(integer_code,translate(aktual_triplet,Alphabete));
	}
	appendValue(trans_reads,integer_code);
	
	

}
