#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1
#include "own_functions.h"



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
// Punkttest f�r g�ltige Eingabeparameter
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
	result = get_translate(triplet,alphabet[0]);	// Glutamins�ure
	SEQAN_ASSERT(result==2);
}
// ung�ltige Eingabeparameter
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
