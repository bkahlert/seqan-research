#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>
#include "read_stats.h"
#include "output.h"
#include "kmer_content.h"

const seqan::Dna5 nucA = 'A';
const seqan::Dna5 nucC = 'C';
const seqan::Dna5 nucG = 'G';
const seqan::Dna5 nucT = 'T';
const seqan::Dna5 nucN = 'N';

const unsigned nucIdxA = ordValue(nucA);
const unsigned nucIdxC = ordValue(nucC);
const unsigned nucIdxG = ordValue(nucG);
const unsigned nucIdxT = ordValue(nucT);
const unsigned nucIdxN = ordValue(nucN);


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
    seqan::CharString qual = "aaaaabcd";
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[1][nucA],0);
   
    stats.collectReadStats(id, "AAAAACGT", qual);
   
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[1][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[2][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[3][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[4][nucA],1);
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucA],0);
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucA],0);
  
    SEQAN_ASSERT_EQ(stats.nucCount[4][nucC],0);
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucC],1);
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucC],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucG],0);
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucG],1);
    SEQAN_ASSERT_EQ(stats.nucCount[7][nucG],0);
    
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucT],0);
    SEQAN_ASSERT_EQ(stats.nucCount[7][nucT],1);
    SEQAN_ASSERT_EQ(stats.nucCount[8][nucT],0);

}


SEQAN_DEFINE_TEST(test_stats_nuc_increments_two_seq_samelength){
    ReadStats stats(9);

    seqan::CharString id = "";
    
    
    seqan::CharString qual = "aaaaaa";


    SEQAN_ASSERT_EQ(stats.nucCount[0][0],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[1][nucIdxA],0);
   
    stats.collectReadStats(id, "CCCAAA", qual);
    stats.collectReadStats(id, "AAACCC", qual);

    SEQAN_ASSERT_GT(5,nucIdxT);
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucIdxT],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucIdxT],0); 
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucIdxA],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucIdxC],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[1][nucIdxA],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[1][nucIdxC],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[2][nucIdxA],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[2][nucIdxC],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[3][nucIdxA],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[3][nucIdxC],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[4][nucIdxA],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[4][nucIdxC],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucIdxA],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucIdxC],1); 

    SEQAN_ASSERT_EQ(stats.nucCount[6][nucIdxA],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucIdxA],0); 

    SEQAN_ASSERT_EQ(stats.nucCount[0][nucIdxG],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucIdxG],0); 
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucIdxN],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucIdxN],0); 
}


SEQAN_DEFINE_TEST(test_stats_nuc_increments_two_seq_getlonger){
    ReadStats stats(9);

    seqan::CharString id = "";


    SEQAN_ASSERT_EQ(stats.nucCount[4][nucIdxA],0);
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucIdxA],0);
    
    stats.collectReadStats(id, "AAAAA", "aaaaa");

    stats.collectReadStats(id, "AAAAAAA", "abcdeaf");
   
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucIdxA],2); 
    SEQAN_ASSERT_EQ(stats.nucCount[1][nucIdxA],2); 
    SEQAN_ASSERT_EQ(stats.nucCount[2][nucIdxA],2); 
    SEQAN_ASSERT_EQ(stats.nucCount[3][nucIdxA],2); 
    SEQAN_ASSERT_EQ(stats.nucCount[4][nucIdxA],2); 
    SEQAN_ASSERT_EQ(stats.nucCount[5][nucIdxA],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucIdxA],1); 
    SEQAN_ASSERT_EQ(stats.nucCount[7][nucIdxA],0); 
    
    SEQAN_ASSERT_EQ(stats.nucCount[0][nucIdxC],0); 
    SEQAN_ASSERT_EQ(stats.nucCount[6][nucIdxC],0); 
    
}

SEQAN_DEFINE_TEST(test_stats_seq_length){
    ReadStats stats(9);

    SEQAN_ASSERT_EQ(0, stats.readLenCount.size()); 
    seqan::CharString id = "";
    
    stats.collectReadStats(id, "AAAAA", "aaaaa");
    stats.collectReadStats(id, "AAAAAAA", "abcdeaf");
    
    SEQAN_ASSERT_EQ(2, stats.readLenCount.size()); 
    SEQAN_ASSERT_EQ(1, stats.readLenCount[7]); 
    SEQAN_ASSERT_EQ(1, stats.readLenCount[5]); 
    
    stats.collectReadStats(id, "GGGGGGG", "abcdeaf");
    
    SEQAN_ASSERT_EQ(2, stats.readLenCount.size()); 
    SEQAN_ASSERT_EQ(2, stats.readLenCount[7]); 
    SEQAN_ASSERT_EQ(1, stats.readLenCount[5]); 


}

SEQAN_DEFINE_TEST(test_output_pos_score){
    ReadStats stats(9);
    seqan::CharString id = "";

    std::string expected = " vgdjagdjkvgdhjsk";

    stats.scoreCount[0]['a'] = 189;
    stats.scoreCount[0]['b'] = 811;
   
    stats.scoreCount[1]['a'] = 500;
    stats.scoreCount[1]['b'] = 500;

    stats.scoreCount[2]['a'] = 1000;

    std::stringstream actual;
    
    TSVWriter w;
    w.writePerPosQual(stats, actual);
    

    SEQAN_ASSERT_EQ(expected, actual.str());

}


SEQAN_DEFINE_TEST(test_collect_quantiles){
    
    ReadStats stats(9);
    
    seqan::String<unsigned> cum_intervals;
    for(unsigned i = 0; i<97;i++){
        seqan::appendValue(cum_intervals,0);
    }

        // start with zeros that the algo has to skip over
        // ... is in the bucket that contains ....
        // mean => 10929
        // lower 10% => 2186
        // lower quartile => 5464
        // upper quartile => 16393
        // upper 10% => 19671

        seqan::appendValue(cum_intervals,35); //a
        seqan::appendValue(cum_intervals,155); //b
        seqan::appendValue(cum_intervals,155); //c
        seqan::appendValue(cum_intervals,155); //d
        seqan::appendValue(cum_intervals,3343); //e lower 10%
        seqan::appendValue(cum_intervals,6455); //f 
        seqan::appendValue(cum_intervals,6458); //g
        seqan::appendValue(cum_intervals,7566); //h lower quart
        seqan::appendValue(cum_intervals,8577); //i
        seqan::appendValue(cum_intervals,9858); //j
        seqan::appendValue(cum_intervals,12377);//k mean
        seqan::appendValue(cum_intervals,16266);//l
        seqan::appendValue(cum_intervals,19377);//m upper quart
        seqan::appendValue(cum_intervals,19397);//n
        seqan::appendValue(cum_intervals,20101);//o upper 10%
        seqan::appendValue(cum_intervals,21001);//p
        seqan::appendValue(cum_intervals,21857);//q


    seqan::String<char> quantiles;

    stats.collectQuantiles(cum_intervals, quantiles, 20);
    SEQAN_ASSERT_EQ(seqan::length(quantiles), 20);

    SEQAN_ASSERT_EQ("eeeffhijjkklllmmmopq", quantiles);
    //mean
    SEQAN_ASSERT_EQ('k', quantiles[9]);

    // lower 10%
    SEQAN_ASSERT_EQ(quantiles[1], 'e');
    
    //upper 10%
    SEQAN_ASSERT_EQ(quantiles[17], 'o');

    // lower quartile
    SEQAN_ASSERT_EQ(quantiles[4], 'f');

    // upper quartile
    SEQAN_ASSERT_EQ(quantiles[14], 'm');


}

SEQAN_DEFINE_TEST(test_collect_gc_per_read_dist){
    
    ReadStats stats(9);
    SEQAN_ASSERT_EQ(10,length(stats.gcPerReadCount));
    
    SEQAN_ASSERT_EQ(0, stats.gcPerReadCount[0]);
    SEQAN_ASSERT_EQ(0, stats.gcPerReadCount[1]);
    SEQAN_ASSERT_EQ(0, stats.gcPerReadCount[9]);

    seqan::CharString id = "";
    stats.collectReadStats(id, "AAAAA", "aaaaa");
   
    // 1 read without any GC
    SEQAN_ASSERT_EQ(1, stats.gcPerReadCount[0]);
    SEQAN_ASSERT_EQ(0, stats.gcPerReadCount[1]);
    SEQAN_ASSERT_EQ(0, stats.gcPerReadCount[9]);

    // add 1 read with 1 C
    stats.collectReadStats(id, "CAAAA", "aaaaa");
    
    SEQAN_ASSERT_EQ(1, stats.gcPerReadCount[0]);
    SEQAN_ASSERT_EQ(1, stats.gcPerReadCount[1]);
    SEQAN_ASSERT_EQ(0, stats.gcPerReadCount[9]);

    // add 1 read with 1 G 
    stats.collectReadStats(id, "GAAAA", "aaaaa");
    
    SEQAN_ASSERT_EQ(1, stats.gcPerReadCount[0]);
    SEQAN_ASSERT_EQ(2, stats.gcPerReadCount[1]);
    SEQAN_ASSERT_EQ(0, stats.gcPerReadCount[9]);

    // add 1 read with all CG 
    stats.collectReadStats(id, "GCGCCCCCC", "aaaaabbbb");
    
    SEQAN_ASSERT_EQ(1, stats.gcPerReadCount[0]);
    SEQAN_ASSERT_EQ(2, stats.gcPerReadCount[1]);
    SEQAN_ASSERT_EQ(1, stats.gcPerReadCount[9]);
} 


SEQAN_DEFINE_TEST(test_collect_per_read_score_dist){
    ReadStats stats(9);
    
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['a']);

    //this will add one mean score of 'j'
    stats.collectReadStats("", "AAAAA", "jjjjj");
   
    SEQAN_ASSERT_EQ(1, stats.perReadQualCount['j']);
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['a']);
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['i']);
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['k']);
    
    //this will add one mean score of 'j'
    stats.collectReadStats("", "AAAAA", "jjjii");
   
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['i']);
    SEQAN_ASSERT_EQ(2, stats.perReadQualCount['j']);
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['a']);

    //this will add one mean score of 'j': fgh -ijk- lmn
    stats.collectReadStats("", "AAAAAA", "fghlmn");
   
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['i']);
    SEQAN_ASSERT_EQ(3, stats.perReadQualCount['j']);
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['a']);
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['f']);
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['h']);
    SEQAN_ASSERT_EQ(0, stats.perReadQualCount['g']);
}

SEQAN_DEFINE_TEST(test_output_nuc_content_dist){
    ReadStats stats(9);
    
    std::stringstream gc_content;
    std::stringstream n_content;
    std::stringstream nuc_dist;
    
    TSVWriter w;
    //w.writeNucContent(stats, gc_content, n_content, nuc_dist);

}

SEQAN_DEFINE_TEST(test_stats_grow_seq_length){
    ReadStats stats(5);

    SEQAN_ASSERT_EQ(length(stats.nucCount),5); 
    
    //this should fit
    stats.collectReadStats("", "AAAAA", "aaaaa");

    //this is too large, it should grow internals to fit: 
    stats.collectReadStats("", "AAAAAAA", "abcdeaf");
        
    SEQAN_ASSERT_EQ(length(stats.nucCount),7); 
   
    
}

SEQAN_DEFINE_TEST(test_kmer_content){
    KmerContent kc(3,3);
   
    kc.CountKmers("ACCATACCA");
   //kc.CountKmers("ATCACGAAT");
         
      
    //SEQAN_ASSERT_EQ(kc.kmerCount["ACC"][0],2); 
    //SEQAN_ASSERT_EQ(kc.kmerCount["CCA"][0],2);
    //SEQAN_ASSERT_EQ(kc.kmerCount["CAT"][0],1);
    //SEQAN_ASSERT_EQ(kc.kmerCount["ATA"][0],1);
    //SEQAN_ASSERT_EQ(kc.kmerCount["TAC"][0],1);

    //SEQAN_ASSERT_EQ(kc.kmerPositions[1]['a'],0);
   
    //stats.collectReadStats(id, seq, qual);
    //SEQAN_ASSERT_EQ(stats.scoreCount[0][0],0); 
    //SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
    //SEQAN_ASSERT_EQ(stats.scoreCount[2]['a'],1);
    //SEQAN_ASSERT_EQ(stats.scoreCount[3]['a'],1);
    //SEQAN_ASSERT_EQ(stats.scoreCount[4]['a'],1);
    //SEQAN_ASSERT_EQ(stats.scoreCount[5]['a'],0);
    //SEQAN_ASSERT_EQ(stats.scoreCount[6]['a'],0);

    
   // SEQAN_ASSERT_EQ(stats.scoreCount[1]['c'],0);
    //SEQAN_ASSERT_EQ(stats.scoreCount[1]['t'],0);
    //SEQAN_ASSERT_EQ(stats.scoreCount[1]['g'],0);
    
    //SEQAN_ASSERT_EQ(stats.scoreCount[1]['a'],1);
    
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

    SEQAN_CALL_TEST(test_stats_seq_length);
    
    SEQAN_CALL_TEST(test_collect_quantiles);
    SEQAN_CALL_TEST(test_collect_gc_per_read_dist);
    SEQAN_CALL_TEST(test_collect_per_read_score_dist);
    SEQAN_CALL_TEST(test_stats_grow_seq_length);
    SEQAN_CALL_TEST(test_kmer_content);

}
SEQAN_END_TESTSUITE


