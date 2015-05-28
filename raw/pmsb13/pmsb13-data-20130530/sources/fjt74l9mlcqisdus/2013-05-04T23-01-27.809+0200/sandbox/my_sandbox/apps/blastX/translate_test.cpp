// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1
#include "own_functions.h"



// TEST HASH FUNKTION --------------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_hash)
{
	// normaler Eingabewert fuer hash funktion
	int result = hash(1,0,3);				
	SEQAN_ASSERT(result>=0 && result<=63);
	// zu hoher Eingabewert in hash funktion
	result = hash(4,7,900);
	SEQAN_ASSERT_EQ(result,-1);
	// zu niedriger Eingabewert in hash funktion
	result = hash(4,-7,-900);
	SEQAN_ASSERT_EQ(result,-1);
	// Eingabewert N
	result = hash(4,0,1);
	SEQAN_ASSERT_EQ(result,-1);
	// cast von normalen char
	result = hash((int)'G',(int)'K',(int)'J');
	SEQAN_ASSERT_EQ(result,-1);
}
// TEST HASH FUNKTION --------------------------------------------------------------------------------------------------


// TEST GET AMINOACID POS FUNKTION -------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_get_amino_acid_pos)
{
	int result;
	// alle Moeglichen Eingabewerte fuer get_amino_acid_pos
	for (int i=0;i<=63;++i){
		result = get_Amino_Acid_Pos(i);
		SEQAN_ASSERT(result>=0 && result<=21);
	}
	// alle Moegliche und nicht Moegliche Eingabewerte
	for (int i=-100;i<=100;++i){
		result = get_Amino_Acid_Pos(i);
		SEQAN_ASSERT(result>=0 && result<=21 || result==-1);
	}
	// Punktprobe
	result = get_Amino_Acid_Pos(7); // Threonin Eingabe
	SEQAN_ASSERT_EQ(result,16);
}
// TEST GET AMINOACID POS FUNKTION -------------------------------------------------------------------------------------


// TEST GET AMINOACID POS UND HASH FUNKTION ----------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_apps_get_amino_acid_pos_and_hash)
{
	// Punkttest hash funktion und get_amino_acid_pos	
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


// TEST GET_TRANSLATE_FROM_CODON ---------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_get_translate_from_codon)
{
	// Punkttest für gültige Eingabeparameter
	int x = 0;
	StrSetSA alphabet = GET_ALPHABET(x,x);

	String<Dna> triplet = "GCT"; 
	int result = get_translate_from_codon(triplet,alphabet[0]);	// Alanin
	SEQAN_ASSERT(result==4);

	triplet = "CTA";
	result = get_translate_from_codon(triplet,alphabet[0]);	// Leucin
	SEQAN_ASSERT(result==4);

	triplet = "GTT";
	result = get_translate_from_codon(triplet,alphabet[0]);	// Valin
	SEQAN_ASSERT(result==4);

	triplet = "TGG";
	result = get_translate_from_codon(triplet,alphabet[0]);	// Tryptophan
	SEQAN_ASSERT(result==5);

	triplet = "GAA";
	result = get_translate_from_codon(triplet,alphabet[0]);	// Glutaminsäure
	SEQAN_ASSERT(result==2);

}
// TEST GET_TRANSLATE_FROM_CODON ---------------------------------------------------------------------------------------



// TEST TRANSLATE ------------------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_translate)
{
	StrSetSA trans_reads;
	int x = 0;
	StrSetSA alphabet = GET_ALPHABET(x,x);
	String<Dna> read = "ACGATGACGATCAGTACGATACAGTAC";
	for (int frame=0;frame<6;++frame){
		int result = translate(trans_reads,read,alphabet[0],frame);
		SEQAN_ASSERT_EQ(result,0);
		SEQAN_ASSERT_EQ(length(trans_reads),frame+1);

	}
	clear(trans_reads);
	read = "";
	for (int frame=0;frame<6;++frame){
		int result = translate(trans_reads,read,alphabet[0],frame);
		SEQAN_ASSERT_EQ(result,0);
		SEQAN_ASSERT_EQ(length(trans_reads),frame+1);

	}
}
// TEST TRANSLATE ------------------------------------------------------------------------------------------------------


// TEST TRANSLATE_READS ------------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_translate_reads)
{
	StrSetSA trans_reads;
	int x = 0;
	StrSetSA alphabet = GET_ALPHABET(x,x);
	StringSet<String<Dna>> read;
	appendValue(read,"ACGATGCAGTCAGTGTCA");
	appendValue(read,"GTGATCGTACGTCAGGTA");
	int result = translate_reads(read,trans_reads,alphabet[0]);
	SEQAN_ASSERT_EQ(length(trans_reads),12);
	
	
	SEQAN_ASSERT_EQ(trans_reads[0],"QEQQCQ");
	SEQAN_ASSERT_EQ(trans_reads[1],"DRQQR");
	SEQAN_ASSERT_EQ(trans_reads[2],"NCCQC");
	/*
	SEQAN_ASSERT_EQ(trans_reads[3],);
	SEQAN_ASSERT_EQ(trans_reads[4],);
	SEQAN_ASSERT_EQ(trans_reads[5],);
	*/

	
}
// TEST TRANSLATE_READS ------------------------------------------------------------------------------------------------


// TEST DURCHLAUF ------------------------------------------------------------------------------------------------------
// Normale eingabe
SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_hash);
	SEQAN_CALL_TEST(test_my_app_get_amino_acid_pos);
	SEQAN_CALL_TEST(test_my_apps_get_amino_acid_pos_and_hash);
	
	SEQAN_CALL_TEST(test_my_app_get_translate_from_codon);
	
	SEQAN_CALL_TEST(test_my_app_translate);
	SEQAN_CALL_TEST(test_my_app_translate_reads);
	//SEQAN_CALL_TEST(test_my_app_translate_database);
	
}
SEQAN_END_TESTSUITE
// TEST DURCHLAUF ------------------------------------------------------------------------------------------------------
