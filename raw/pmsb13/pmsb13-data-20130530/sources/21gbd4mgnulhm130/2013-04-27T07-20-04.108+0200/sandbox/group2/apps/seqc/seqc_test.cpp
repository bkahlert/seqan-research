#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>
#include "read_stats.h"

const unsigned nucA = (seqan::Dna5)'A';
const unsigned nucC = (seqan::Dna5)'C';
const unsigned nucG = (seqan::Dna5)'G';
const unsigned nucT = (seqan::Dna5)'T';
const unsigned nucN = (seqan::Dna5)'N';

SEQAN_DEFINE_TEST(test_stats_inits_nuc_count_to_zero){
    
    ReadStats stats(9);

    SEQAN_ASSERT_EQ(stats.nucCount[0][nucA],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][nucA],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucC],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][nucC],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucG],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][nucG],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucT],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][nucT],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucT],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[8][nucT],0);
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


SEQAN_DEFINE_TEST(test_stats_score_increments_two_seq_samelength){
    ReadStats stats(9);

    seqan::CharString id = "";
    seqan::Dna5String seq;
    
    
    seq = "ACGTA";
    seqan::CharString qual1 = "aaaaa";
    seqan::CharString qual2 = "abcde";


    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],0);
   
    stats.collectReadStats(id, seq, qual1);

    stats.collectReadStats(id, seq, qual2);
    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 

    SEQAN_ASSERT_EQ(stats.scoreCount[0]['a'],2);
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
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
    
}


SEQAN_DEFINE_TEST(test_stats_score_increments_two_seq_getlonger){
    ReadStats stats(9);

    seqan::CharString id = "";


    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],0);
   
    stats.collectReadStats(id, "AAAAA", "aaaaa");

    stats.collectReadStats(id, "AAAAAAA", "abcdeaf");
    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[0]['a'],2);
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
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
    
    
    SEQAN_ASSERT_EQ(stats.scoreCount[3]['e'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[4]['e'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[5]['e'],0);
    
    SEQAN_ASSERT_EQ(stats.scoreCount[4]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[5]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[6]['a'],0);

    SEQAN_ASSERT_EQ(stats.scoreCount[5]['f'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[6]['f'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[7]['f'],0);
    
}

SEQAN_DEFINE_TEST(test_stats_nuc_increments_single_seq){
    ReadStats stats(9);

    seqan::CharString id = "";
    seqan::CharString qual = "aaaaa";
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[1][nucA],0);
   
    stats.collectReadStats(id, "AAAAA", qual);
   
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[1][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[2][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[3][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[4][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucA],0);
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucA],0);

}


SEQAN_DEFINE_TEST(test_stats_nuc_increments_two_seq_samelength){
    ReadStats stats(9);

    seqan::CharString id = "";
    seqan::Dna5String seq;
    
    
    seq = "ACGTA";
    seqan::CharString qual1 = "aaaaa";
    seqan::CharString qual2 = "abcde";


    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],0);
   
    stats.collectReadStats(id, seq, qual1);

    stats.collectReadStats(id, seq, qual2);
    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 

    SEQAN_ASSERT_EQ(stats.scoreCount[0]['a'],2);
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
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
    
}


SEQAN_DEFINE_TEST(test_stats_nuc_increments_two_seq_getlonger){
    ReadStats stats(9);

    seqan::CharString id = "";


    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],0);
   
    stats.collectReadStats(id, "AAAAA", "aaaaa");

    stats.collectReadStats(id, "AAAAAAA", "abcdeaf");
    SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.scoreCount[0]['a'],2);
    SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
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
    
    
    SEQAN_ASSERT_EQ(stats.scoreCount[3]['e'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[4]['e'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[5]['e'],0);
    
    SEQAN_ASSERT_EQ(stats.scoreCount[4]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[5]['a'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[6]['a'],0);

    SEQAN_ASSERT_EQ(stats.scoreCount[5]['f'],0);
    SEQAN_ASSERT_EQ(stats.scoreCount[6]['f'],1);
    SEQAN_ASSERT_EQ(stats.scoreCount[7]['f'],0);
    
}






SEQAN_BEGIN_TESTSUITE(test_seqc_stats)
{
    SEQAN_CALL_TEST(test_stats_inits_score_count_to_zero);
    SEQAN_CALL_TEST(test_stats_inits_nuc_count_to_zero);
    
    SEQAN_CALL_TEST(test_stats_score_increments_single_seq);
    SEQAN_CALL_TEST(test_stats_score_increments_two_seq_samelength);
    SEQAN_CALL_TEST(test_stats_score_increments_two_seq_getlonger);

    SEQAN_CALL_TEST(test_stats_nuc_increments_single_seq);
    SEQAN_CALL_TEST(test_stats_nuc_increments_two_seq_samelength);
    SEQAN_CALL_TEST(test_stats_nuc_increments_two_seq_getlonger);

}
SEQAN_END_TESTSUITE


