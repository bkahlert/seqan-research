#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include "output.h"
#include "read_stats.h"

/*
 * Text constants for output
 */

// output filenames
static const char * SEQC_READLENGTH_FILE        = "readlength.tsv";
static const char * SEQC_POS_QUAL_FILE          = "per_pos_qual.tsv";
static const char * SEQC_POS_GC_CONTENT_FILE    = "pos_gc_content.tsv";
static const char * SEQC_POS_N_CONTENT_FILE     = "pos_n_content.tsv";
static const char * SEQC_POS_NUC_DIST_FILE      = "pos_nuc_dist.tsv";
static const char * SEQC_SEQ_GC_CONTENT_FILE    = "seq_gc_content.tsv";
static const char * SEQC_SEQ_QUAL_DIST_FILE     = "seq_qual.tsv";
static const char * SEQC_SUMMARY_FILE           = "summary.tsv";
 
// textual conventions
static const char * SEQC_END_MODULE_LN = ">>END_MODULE"; 

/**
 * wrapper: actually writes statistics info to files. 
 */
int TSVWriter::writeAll(ReadStats & stats, std::string const outdir){
    if (mkdir("testout",0744) != 0) {
         std::cerr <<  "Error creating output directory: " << strerror( errno ) << std::endl;
        
        exit(-1);
    } 

    std::ofstream ofStream;
    
    {// ReadLength Distribution
    ofStream.open((outdir + "/" + SEQC_READLENGTH_FILE).c_str());
    std::stringstream ss;
    writeReadLength(stats, ss);    
    ofStream << ss.str();
    ofStream.close();
    }


    {// Per Position-in-Read Quality Statistics
    ofStream.open((outdir + "/" + SEQC_POS_QUAL_FILE).c_str());
    std::stringstream ss;
    writePerPosQual(stats, ss);    
    ofStream << ss.str();
    ofStream.close();
    }

    {// Per Position-in-Read Nucleotide content Statistics
    std::stringstream gc_content;
    std::stringstream n_content;
    std::stringstream nuc_dist;
    writeNucContent(stats, gc_content, n_content, nuc_dist);    
    
    ofStream.open((outdir + "/" + SEQC_POS_GC_CONTENT_FILE).c_str());
    ofStream << gc_content.str();
    ofStream.close();

    std::ofstream ofStream2;
    ofStream2.open((outdir + "/" + SEQC_POS_N_CONTENT_FILE).c_str());
    ofStream2 << n_content.str();
    ofStream2.close();
    
    std::ofstream ofStream3;
    ofStream3.open((outdir + "/" + SEQC_POS_NUC_DIST_FILE).c_str());
    ofStream3 << nuc_dist.str();
    ofStream3.close();
    }


}

/**
 * create the content that goes into the sequence length distribution statistics
 */
void TSVWriter::writeReadLength(ReadStats & stats, std::stringstream & ss)
{    
    // write out stats
    ss << ">>Sequence Length Distribution\tpass" << std::endl;
    ss << "#Length\tCount" << std::endl;
    std::map<int,int>::iterator iter; 
  
    for(iter=stats.readLenCount.begin();iter != stats.readLenCount.end(); iter++){
        ss << iter->first << "\t" << iter->second << std::endl;
    }
    
    ss << SEQC_END_MODULE_LN << std::endl; 

};

unsigned sum(seqan::String<unsigned >  counts ){
    unsigned sum = 0; 
    for(unsigned c = 0; c < length(counts);++c){
        sum += counts[c];
    }
    return sum;

}

void TSVWriter::writeNucContent(ReadStats & stats, std::stringstream & gc_content, 
                                std::stringstream & n_content, std::stringstream & nuc_dist )
{
    gc_content << ">>Per base GC content\tpass" << std::endl;
    gc_content << "#Base\t%GC" << std::endl;
    
    n_content << ">>Per base N content\tpass" << std::endl;
    n_content << "#Base\t%GC" << std::endl;

    nuc_dist << ">>Per base sequence content\tpass" << std::endl;
    nuc_dist << "#Base\tG\tA\tT\tC\tN" << std::endl;
    
    //per position in read:
    for(unsigned pos = 0; pos < length(stats.nucCount); ++pos)
    {

        unsigned sum_a = 0;
        unsigned sum_t = 0;
        unsigned sum_c = 0;
        unsigned sum_g = 0;
        unsigned sum_n = 0;
        
        sum_c += sum(stats.nucCount[pos][(seqan::Dna5)'C']);
        sum_g += sum(stats.nucCount[pos][(seqan::Dna5)'G']);
        sum_a += sum(stats.nucCount[pos][(seqan::Dna5)'A']);
        sum_t += sum(stats.nucCount[pos][(seqan::Dna5)'T']);
        sum_n += sum(stats.nucCount[pos][(seqan::Dna5)'N']);

        unsigned sum_all = sum_c + sum_g + sum_a + sum_t + sum_n;

        if(sum_all > 0)
        {
            gc_content << (pos+1) << "\t"; 
            gc_content << (float)(sum_g+sum_c)/sum_all  << std::endl;

            n_content << (pos+1) << "\t";
            n_content << (float)sum_n/sum_all << std::endl;

            nuc_dist << (pos+1) << "\t";
            nuc_dist 
                << (float)sum_g/sum_all  << "\t"
                << (float)sum_a/sum_all  << "\t"
                << (float)sum_t/sum_all  << "\t"
                << (float)sum_c/sum_all  << "\t"
                << (float)sum_n/sum_all 
                << std::endl;
            
        }

    }
};

void TSVWriter::writePerPosQual(ReadStats & stats, std::stringstream & ss)
{
    unsigned score_offset = 33;
    ss << ">>Per base sequence quality\tfail" << std::endl;
    ss << "Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile" 
        << std::endl;
  
            //per position in read:
            for(unsigned pos = 0; pos < length(stats.scoreCount); ++pos){
                // we have to know how many there are counted in total
                // so we can pick specific quantiles later
                unsigned total_scores = 0;
                
                //we'll cumulate counts here so we have easier access later  
                seqan::String<unsigned> intervals;
               
                //we need to add all scores to calc mean
                unsigned score_sum = 0;
              
                //we'll collect 20 quantiles here
                seqan::String<char> quantiles;

                for(unsigned score = 0; score < length(stats.scoreCount[0]); ++score){
                    total_scores += stats.scoreCount[pos][score];
                    seqan::append(intervals,total_scores);
                    score_sum += stats.scoreCount[pos][score]*score;
                }
                if(total_scores > 0){                
                    stats.collectQuantiles(intervals, quantiles, 20);
                    
                    //current base pos and mean
                    ss << (pos+1) << "\t" << ((double)score_sum/total_scores)-score_offset ;                
                    ss << "\t" <<  quantiles[9] - score_offset;  // median
                    ss << "\t" <<  quantiles[4] - score_offset;  // lower quartile
                    ss << "\t" <<  quantiles[14] - score_offset;  // upper quartile
                    ss << "\t" <<  quantiles[1] - score_offset;  // 10th percentile
                    ss << "\t" <<  quantiles[17] - score_offset;  // 90th percentile
                    ss << std::endl;
                }
        }
    


    //"
};
