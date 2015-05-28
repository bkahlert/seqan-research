// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1
#include "own_functions.h"

/*
// TEST GET_ALPHABET_FORCE ---------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_get_alphabet_force)
{
	StringSet<String<AminoAcid> > result;	
	String<AminoAcid> temp;
	append(temp, "CDQNRQNCDCCDRCCQQQQCE");
	appendValue(result, temp);
	StringSet<String<AminoAcid> > alphabet = GET_ALPHABET_FORCE();
	//SEQAN_ASSERT_EQ(result, alphabet);
}
// TEST GET_ALPHABET_FORCE ---------------------------------------------------------------------------------------------


// TEST GET_ALPHABET ---------------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_get_alphabet)
{

}
// TEST GET_ALPHABET ---------------------------------------------------------------------------------------------------
*/

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
	StrSetSA alphabet = GET_ALPHABET_FORCE();

	String<Dna> triplet = "GCT"; 
	int result = get_translate_from_codon(triplet,alphabet[0],1);	// Alanin
	SEQAN_ASSERT(result==4);

	triplet = "CTA";
	result = get_translate_from_codon(triplet,alphabet[0],1);	// Leucin
	SEQAN_ASSERT(result==4);

	triplet = "GTT";
	result = get_translate_from_codon(triplet,alphabet[0],1);	// Valin
	SEQAN_ASSERT(result==4);

	triplet = "TGG";
	result = get_translate_from_codon(triplet,alphabet[0],1);	// Tryptophan
	SEQAN_ASSERT(result==5);

	triplet = "GAA";
	result = get_translate_from_codon(triplet,alphabet[0],1);	// Glutaminsäure
	SEQAN_ASSERT(result==2);

	triplet = "GCT";
	result = get_translate_from_codon(triplet,alphabet[0],0);	// Alanin
	SEQAN_ASSERT(result==0);

	triplet = "CTA";
	result = get_translate_from_codon(triplet,alphabet[0],0);	// Leucin
	SEQAN_ASSERT(result==10);

	triplet = "GTT";
	result = get_translate_from_codon(triplet,alphabet[0],0);	// Valin
	SEQAN_ASSERT(result==19);

	triplet = "TGG";
	result = get_translate_from_codon(triplet,alphabet[0],0);	// Tryptophan
	SEQAN_ASSERT(result==17);

	triplet = "GAA";
	result = get_translate_from_codon(triplet,alphabet[0],0);	// Glutaminsäure
	SEQAN_ASSERT(result==6);

}
// TEST GET_TRANSLATE_FROM_CODON ---------------------------------------------------------------------------------------



// TEST TRANSLATE ------------------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_translate)
{
	StrSetSA trans_reads;
	StrSetSA alphabet = GET_ALPHABET_FORCE();
	String<Dna> read = "ACGATGACGATCAGTACGATACAGTAC";
	for (int frame=0;frame<6;++frame){
		int result = translate(trans_reads,read,alphabet[0],frame,1);
		SEQAN_ASSERT_EQ(result,0);
		SEQAN_ASSERT_EQ(length(trans_reads),frame+1);

	}
	SEQAN_ASSERT_EQ(trans_reads[0],"QRQCQQCQQ");
	SEQAN_ASSERT_EQ(trans_reads[1],"DEDQCDQQ");
	SEQAN_ASSERT_EQ(trans_reads[2],"NNNQQNQC");
	SEQAN_ASSERT_EQ(trans_reads[3],"CCQDQNDDD");
	SEQAN_ASSERT_EQ(trans_reads[4],"QRCCCCCC");
	SEQAN_ASSERT_EQ(trans_reads[5],"QCQQEQQQ");
	clear(trans_reads);	
	for (int frame=0;frame<6;++frame){
		int result = translate(trans_reads,read,alphabet[0],frame,0);
		SEQAN_ASSERT_EQ(result,0);
		SEQAN_ASSERT_EQ(length(trans_reads),frame+1);
	}
	SEQAN_ASSERT_EQ(trans_reads[0],"TMTISTIQY");
	SEQAN_ASSERT_EQ(trans_reads[1],"RBRSVRYS");
	SEQAN_ASSERT_EQ(trans_reads[2],"DDDQYDTV");
	SEQAN_ASSERT_EQ(trans_reads[3],"VLYRTDRHR");
	SEQAN_ASSERT_EQ(trans_reads[4],"YCIVLIVI");
	SEQAN_ASSERT_EQ(trans_reads[5],"TVSYBSSS");
	
}
// TEST TRANSLATE ------------------------------------------------------------------------------------------------------


// TEST TRANSLATE_READS ------------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_translate_reads)
{
	StrSetSA trans_reads;
	StrSetSA alphabet = GET_ALPHABET_FORCE();
	StringSet<String<Dna>> read;
	appendValue(read,"ACGATGCAGTCAGTGTCA");
	appendValue(read,"GTGATCGTACGTCAGGTA");
	int result = translate_reads(read,trans_reads,alphabet[0],1);
	SEQAN_ASSERT_EQ(length(trans_reads),12);
	
	
	SEQAN_ASSERT_EQ(trans_reads[0],"QRQQCQ");
	SEQAN_ASSERT_EQ(trans_reads[1],"DRQQR");
	SEQAN_ASSERT_EQ(trans_reads[2],"NCCQC");
	
	SEQAN_ASSERT_EQ(trans_reads[3],"EDECDD");
	SEQAN_ASSERT_EQ(trans_reads[4],"NQNRC");
	SEQAN_ASSERT_EQ(trans_reads[5],"QCQCQ");

	
	SEQAN_ASSERT_EQ(trans_reads[6],"CCCDQC");
	SEQAN_ASSERT_EQ(trans_reads[7],"EQQCD");
	SEQAN_ASSERT_EQ(trans_reads[8],"NDQQC");

	SEQAN_ASSERT_EQ(trans_reads[9],"QCQQND");
	SEQAN_ASSERT_EQ(trans_reads[10],"QEDQC");
	SEQAN_ASSERT_EQ(trans_reads[11],"CNCDQ");
}
// TEST TRANSLATE_READS ------------------------------------------------------------------------------------------------

/*
// TEST TRANSLATE_DATABASE ---------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_translate_database)
{
	StrSetSA trans_proteine;
	StrSetSA alphabet = GET_ALPHABET_FORCE();
	StrSetSA proteine;
	appendValue(proteine,"QNEVGNATILNVWPF");
	appendValue(proteine,"ARNDCQEGHILKMFPSTWYV");
	int result = translate_database(trans_proteine,proteine,alphabet[0]);
	SEQAN_ASSERT_EQ(length(trans_proteine),2);
	SEQAN_ASSERT_EQ(trans_proteine[0],"QQNCCQCQCCQCQCC");
	SEQAN_ASSERT_EQ(trans_proteine[1],"CDQNRQNCDCCDRCCQQQQC");
}
// TEST TRANSLATE_DATABASE ---------------------------------------------------------------------------------------------

// TEST APPEND_TO_MATCH_FOUND ------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_append_to_match_found)
{
	Match_found test;
	Pair<unsigned int,unsigned int> x (1,10);
	Pair<unsigned int,unsigned int> y (1,13);
	Pair<unsigned int,unsigned int> z (13,0);
	int seed = 3;
	append_to_match_found(test,x,y,z,seed,10);
	SEQAN_ASSERT_EQ(test.position_read,1);
	SEQAN_ASSERT_EQ(test.begin_read,5);
	SEQAN_ASSERT_EQ(test.end_read,8);
	SEQAN_ASSERT_EQ(test.position_protein,1);
	SEQAN_ASSERT_EQ(test.begin_protein,5);
	SEQAN_ASSERT_EQ(test.end_protein,8);

}
// TEST APPEND_TO_MATCH_FOUND ------------------------------------------------------------------------------------------

// TEST APPEND_TO_FIND_MATCHES -----------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_find_matches)
{
	StrSetSA trans_proteine;
	StrSetSA trans_reads;
	appendValue (trans_reads,"QCQ");
	appendValue(trans_proteine,"CRCQQCQQ");
	Match_found test0;
	int seed = 3;
	int distance = 0;
	int result = find_matches(trans_proteine,trans_reads,seed,distance,test0);
	SEQAN_ASSERT_EQ(result,0);
	SEQAN_ASSERT_EQ(test0.position_read,0);
	SEQAN_ASSERT_EQ(test0.begin_read,0);
	SEQAN_ASSERT_EQ(test0.end_read,3);
	SEQAN_ASSERT_EQ(test0.position_protein,0);
	SEQAN_ASSERT_EQ(test0.begin_protein,4);
	SEQAN_ASSERT_EQ(test0.end_protein,7);

	appendValue (trans_reads,"CQQ");
	appendValue(trans_proteine,"QCQQCRCQ");
	Match_found test;
	result = find_matches(trans_proteine,trans_reads,seed,distance,test);

}
// TEST APPEND_TO_FIND_MATCHES -----------------------------------------------------------------------------------------

// TEST APPEND_TO_FIND_MATCHES_FOR_ALL ---------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_find_matches_for_all)
{

}
// TEST APPEND_TO_FIND_MATCHES_FOR_ALL ---------------------------------------------------------------------------------

*/

// TEST DURCHLAUF ------------------------------------------------------------------------------------------------------
// Normale eingabe
SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
   // SEQAN_CALL_TEST(test_my_app_get_alphabet_force);
//	SEQAN_CALL_TEST(test_my_app_get_alphabet);
	
	SEQAN_CALL_TEST(test_my_app_hash);
	SEQAN_CALL_TEST(test_my_app_get_amino_acid_pos);
	SEQAN_CALL_TEST(test_my_apps_get_amino_acid_pos_and_hash);
	SEQAN_CALL_TEST(test_my_app_get_translate_from_codon);
	SEQAN_CALL_TEST(test_my_app_translate);
	SEQAN_CALL_TEST(test_my_app_translate_reads);
//	SEQAN_CALL_TEST(test_my_app_translate_database);
	
//	SEQAN_CALL_TEST(test_my_app_append_to_match_found);
//	SEQAN_CALL_TEST(test_my_app_find_matches);
//	SEQAN_CALL_TEST(test_my_app_find_matches_for_all);

	
	

	
	
}
SEQAN_END_TESTSUITE
// TEST DURCHLAUF ------------------------------------------------------------------------------------------------------
