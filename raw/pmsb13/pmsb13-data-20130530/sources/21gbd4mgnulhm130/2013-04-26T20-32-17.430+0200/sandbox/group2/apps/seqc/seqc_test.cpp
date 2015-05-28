#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>
#include "read_stats.h"


SEQAN_DEFINE_TEST(test_stats_inits_nuc_count_to_zero){
    ReadStats stats(9);
    unsigned a = (seqan::Dna5)'A';
    unsigned c = (seqan::Dna5)'C';
    unsigned g = (seqan::Dna5)'G';
    unsigned t = (seqan::Dna5)'T';
    unsigned n = (seqan::Dna5)'N';

    SEQAN_ASSERT_EQ(stats.nucCount[0][a],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][a],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][c],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][c],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][g],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][g],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][t],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][t],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][n],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][n],0);
} 

SEQAN_DEFINE_TEST(test_stats_inits_score_count_to_zero){
    ReadStats stats(9);
    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[8][8],0);
    
}


SEQAN_DEFINE_TEST(test_stats_score_increments_single_seq){
    ReadStats stats(9);

    seqan::CharString id = "";
    seqan::Dna5String seq;
    seqan::CharString qual;
    
    
    seq = "ACGTA";
    qual = "aaaaa";


    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],0);
   
    stats.collectReadStats(id, seq, qual);
    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[2]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[3]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[4]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[5]['a'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[6]['a'],0);

    
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['c'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['t'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['g'],0);
    
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
    
}


SEQAN_DEFINE_TEST(test_stats_score_increments_two_seq){
    ReadStats stats(9);

    seqan::CharString id = "";
    seqan::Dna5String seq;
    seqan::CharString qual;
    
    
    seq = "ACGTA";
    qual = "aaaaa";


    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],0);
   
    stats.collectReadStats(id, seq, qual);

    qual = "abcd";
    stats.collectReadStats(id, seq, qual);
    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
/*
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],2);
    SEQAN_ASSERT_EQ(stats.scoreCount[2]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[3]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[4]['a'],1);

    SEQAN_ASSERT_EQ(stats.scoreCount[0]['b'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['b'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[2]['b'],0);

    SEQAN_ASSERT_EQ(stats.scoreCount[1]['c'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[2]['c'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[3]['c'],0);
    
    SEQAN_ASSERT_EQ(stats.scoreCount[2]['d'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[3]['d'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[4]['d'],0);
    
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
  */  
}




SEQAN_BEGIN_TESTSUITE(test_seqc)
{
    SEQAN_CALL_TEST(test_stats_inits_score_count_to_zero);
    SEQAN_CALL_TEST(test_stats_inits_nuc_count_to_zero);
    SEQAN_CALL_TEST(test_stats_score_increments_single_seq);
    SEQAN_CALL_TEST(test_stats_score_increments_two_seq);

}

SEQAN_END_TESTSUITE
