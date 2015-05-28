#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1
#include "own_functions.h"



// GETESTETE FUNKTIONEN ------------------------------------------------------------------------------------------------
int hash(int x,int y,int z){
	int erg = (x*16)+(y*4)+z;
	if (erg>=0 && erg<=63) return erg;
	else {
		cerr << "invalid input in hash-function" <<endl;
		return -1;
	}
}
// GETESTETE FUNKTIONEN ------------------------------------------------------------------------------------------------
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
// GETESTETE FUNKTIONEN ------------------------------------------------------------------------------------------------
int get_translate(String<Dna5> & aktual_triplet,String<AminoAcid> & Alphabete){
	for (int i=0;i<length(aktual_triplet);++i){
		if (aktual_triplet[i]<0 && aktual_triplet[i]>3) return -1;
	}
	int hash_value = hash((int) aktual_triplet[0],(int) aktual_triplet[1],(int) aktual_triplet[2]);
	if (hash_value!=-1){
		int Amino_pos = get_Amino_Acid_Pos(hash_value);
		if (Amino_pos!=-1){
			int gruppe = Alphabete[Amino_pos];
			return gruppe;
		}
		else return -1;
	}
	else return -1;
}
// GETESTETE FUNKTIONEN ------------------------------------------------------------------------------------------------
StrSetSA GET_ALPHABET(int & numb_alp,int & size_alp){
	StrSetSA Alphabete;
	String<AminoAcid>alp1; 
	append(alp1,4); // Alanin			0
	append(alp1,3); // Arginin			1
	append(alp1,5); // Asparagin		2
	append(alp1,2); // Asparaginsäure	3
	append(alp1,1); // Cystein			4
	append(alp1,5);	// Glutamin			5
	append(alp1,2); // Glutaminsäure	6
	append(alp1,4); // Glycin			7
	append(alp1,3); // Histidin			8
	append(alp1,4); // Isoleucin		9
	append(alp1,4); // Leucin			10
	append(alp1,3); // Lysin			11
	append(alp1,1); // Methionin		12
	append(alp1,4); // Phenylalanin		13
	append(alp1,4);	// Prolin			14
	append(alp1,5); // Serin			15
	append(alp1,5); // Threonin			16
	append(alp1,5); // Tryptophan		17
	append(alp1,5); // Tyrosin			18
	append(alp1,4); // Valin			19
	append(alp1,6);	// STOP				20
	appendValue(Alphabete,alp1);
	return (Alphabete);
}

// GETESTETE FUNKTIONEN ------------------------------------------------------------------------------------------------

int translate_reads(StrSetSA & trans_reads,String<Dna5> & read,String<AminoAcid> & Alphabete,int frame){
	String<Dna5>revComplement;
	// frame = 0,1 oder 2 steht fuer forward 
	// frame = 3,4,oder 5 steht fuer reverse, wenn also frame/3==1 dann wird das
	// complement gebildet und alle frames dafuer durchgegangen 
	if ((frame/3)==1){
		revComplement=read;
		reverseComplement(revComplement);
		frame -= 3;
	}
	String<AminoAcid>integer_code;
	// Schleife durch Reads mit Frame Verschiebung
	for (int triplet=0+frame;triplet<length(read) && triplet+3<=length(read);triplet+=3){
		String<Dna5> aktual_triplet;
		if ((frame/3)==0) aktual_triplet = infix(read,triplet,triplet+3);
		else aktual_triplet = infix(revComplement,triplet,triplet+3);
		int translate = get_translate(aktual_triplet,Alphabete);
		if (translate!=-1) appendValue(integer_code,translate);
		else return 1;
	}
	appendValue(trans_reads,integer_code);
	return 0;
}
// GETESTETE FUNKTIONEN ------------------------------------------------------------------------------------------------






// TEST HASH FUNKTION --------------------------------------------------------------------------------------------------
// normaler Eingabewert fuer hash funktion
SEQAN_DEFINE_TEST(test_my_app_hash_normal)
{
	int result = hash(1,0,3);
	SEQAN_ASSERT(result>=0 && result<=63);
}
// zu hoher Eingabewert in hash funktion
SEQAN_DEFINE_TEST(test_my_app_hash_zu_hoch)
{
	int result = hash(4,7,900);
	SEQAN_ASSERT_EQ(result,-1);
}
// zu niedriger Eingabewert in hash funktion
SEQAN_DEFINE_TEST(test_my_app_hash_zu_niedrig)
{
	int result = hash(4,-7,-900);
	SEQAN_ASSERT_EQ(result,-1);
}
// knapp falsche Eingabewerte
SEQAN_DEFINE_TEST(test_my_app_hash_knapp_zu_hoch)
{
	int result = hash(3,3,4);
	SEQAN_ASSERT_EQ(result,-1);
}
// casten von anderen Buchstaben
SEQAN_DEFINE_TEST(test_my_app_hash_cast)
{
	cout <<(int)'G' <<"\t"<<(int)'K' <<"\t"<<(int)'J'<<endl;
	int result = hash((int)'G',(int)'K',(int)'J');
	SEQAN_ASSERT_EQ(result,-1);
}
// TEST HASH FUNKTION --------------------------------------------------------------------------------------------------

// TEST GET AMINOACID POS FUNKTION -------------------------------------------------------------------------------------
// alle Moeglichen Eingabewerte fuer get_amino_acid_pos
SEQAN_DEFINE_TEST(test_my_app_get_amino_acid_pos)
{
	for (int i=0;i<=63;++i){
		int result = get_Amino_Acid_Pos(i);
		SEQAN_ASSERT(result>=0 && result<=21);
	}
}
// alle Moegliche und nicht Moegliche Eingabewerte
SEQAN_DEFINE_TEST(test_my_app_get_amino_acid_pos_outside)
{
	for (int i=-100;i<=100;++i){
		int result = get_Amino_Acid_Pos(i);
		SEQAN_ASSERT(result>=0 && result<=21 || result==-1);
	}
}
// TEST GET AMINOACID POS FUNKTION -------------------------------------------------------------------------------------

// TEST GET AMINOACID POS UND HASH FUNKTION ----------------------------------------------------------------------------
// Punkttest hash funktion und get_amino_acid_pos
SEQAN_DEFINE_TEST(test_my_apps_get_amino_acid_pos_and_hash)
{
	int hash_value = hash(0,0,0);	// Lysin
	int result = get_Amino_Acid_Pos(hash_value);	
	SEQAN_ASSERT(result==11);

	hash_value = hash(3,3,3);	// Phenylalanin
	result = get_Amino_Acid_Pos(hash_value);	
	SEQAN_ASSERT(result==13);

	hash_value = hash(1,0,3);	// Histidin
	result = get_Amino_Acid_Pos(hash_value);	
	SEQAN_ASSERT(result==8);

	hash_value = hash(0,2,3);	// Serin
	result = get_Amino_Acid_Pos(hash_value);	
	SEQAN_ASSERT(result==15);

	hash_value = hash(3,2,0);	// stop
	result = get_Amino_Acid_Pos(hash_value);	
	SEQAN_ASSERT(result==20);
	
}
// TEST GET AMINOACID POS UND HASH FUNKTION ----------------------------------------------------------------------------

// TEST GET_TRANSLATE --------------------------------------------------------------------------------------------------
// Punkttest für gültige Eingabeparameter
SEQAN_DEFINE_TEST(test_my_app_get_translate)
{
	int x = 0;
	StrSetSA alphabet = GET_ALPHABET(x,x);

	String<Dna5> triplet = "GCT"; 
	int result = get_translate(triplet,alphabet[0]);	// Alanin
	SEQAN_ASSERT(result==4);

	triplet = "CTA";
	result = get_translate(triplet,alphabet[0]);	// Leucin
	SEQAN_ASSERT(result==4);

	triplet = "GTT";
	result = get_translate(triplet,alphabet[0]);	// Valin
	SEQAN_ASSERT(result==4);

	triplet = "TGG";
	result = get_translate(triplet,alphabet[0]);	// Tryptophan
	SEQAN_ASSERT(result==5);

	triplet = "GAA";
	result = get_translate(triplet,alphabet[0]);	// Glutaminsäure
	SEQAN_ASSERT(result==2);
}
// ungültige Eingabeparameter
SEQAN_DEFINE_TEST(test_my_app_get_translate_invalid)
{
	int x = 0;
	StrSetSA alphabet = GET_ALPHABET(x,x);

	String<Dna5> triplet = "JHI"; 
	int result = get_translate(triplet,alphabet[0]);
	SEQAN_ASSERT_EQ(result,-1);

	triplet = "OMK";
	result = get_translate(triplet,alphabet[0]);	
	SEQAN_ASSERT_EQ(result,-1);

	triplet = "NTT";
	result = get_translate(triplet,alphabet[0]);	
	SEQAN_ASSERT_EQ(result,-1);
}
// TEST GET_TRANSLATE --------------------------------------------------------------------------------------------------

// TEST TRANSLATE_READS ------------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_translate_reads)
{
	StrSetSA trans_reads;
	int x = 0;
	StrSetSA alphabet = GET_ALPHABET(x,x);
	String<Dna5> read = "ACGATGACGATCAGTACGATACAGTAC";
	for (int frame=0;frame<6;++frame){
		int result = translate_reads(trans_reads,read,alphabet[0],frame);
		SEQAN_ASSERT_EQ(result,0);
		SEQAN_ASSERT_EQ(length(trans_reads),frame+1);

	}
}
//
SEQAN_DEFINE_TEST(test_my_app_translate_reads_leer)
{
	StrSetSA trans_reads;
	int x = 0;
	StrSetSA alphabet = GET_ALPHABET(x,x);
	String<Dna5> read = "";
	for (int frame=0;frame<6;++frame){
		int result = translate_reads(trans_reads,read,alphabet[0],frame);
		SEQAN_ASSERT_EQ(result,0);
		SEQAN_ASSERT_EQ(length(trans_reads),frame+1);

	}
}//
SEQAN_DEFINE_TEST(test_my_app_translate_reads_falsch_read)
{
	StrSetSA trans_reads;
	int x = 0;
	StrSetSA alphabet = GET_ALPHABET(x,x);
	String<Dna5> read = "ACGNACGANCAGNAGCTGAHSJJOK";
	for (int frame=0;frame<6;++frame){
		int result = translate_reads(trans_reads,read,alphabet[0],frame);
		SEQAN_ASSERT_EQ(result,1);
	}
}
// TEST TRANSLATE_READS ------------------------------------------------------------------------------------------------




// TEST DURCHLAUF ------------------------------------------------------------------------------------------------------
// Normale eingabe
SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_hash_normal);
	SEQAN_CALL_TEST(test_my_app_hash_zu_hoch);
	SEQAN_CALL_TEST(test_my_app_hash_zu_niedrig);
	SEQAN_CALL_TEST(test_my_app_hash_knapp_zu_hoch);
	SEQAN_CALL_TEST(test_my_app_hash_cast);
	
	SEQAN_CALL_TEST(test_my_app_get_amino_acid_pos);
	SEQAN_CALL_TEST(test_my_app_get_amino_acid_pos_outside);
	
	SEQAN_CALL_TEST(test_my_apps_get_amino_acid_pos_and_hash);
	
	SEQAN_CALL_TEST(test_my_app_get_translate);
	SEQAN_CALL_TEST(test_my_app_get_translate_invalid);
	
	SEQAN_CALL_TEST(test_my_app_translate_reads);
	SEQAN_CALL_TEST(test_my_app_translate_reads_leer);
	SEQAN_CALL_TEST(test_my_app_translate_reads_falsch_read);
}
SEQAN_END_TESTSUITE
// TEST DURCHLAUF ------------------------------------------------------------------------------------------------------
