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

void TSVWriter::writePerPosQual(ReadStats & stats, std::stringstream & ss)
{
    ss << ">>Per base sequence quality\tfail" << std::endl;
    ss << "Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile" 
        << std::endl;
   
        for(unsigned i = 0; i < length(stats.scoreCount); ++i)
            ss <<  i <<": "  << stats.scoreCount[i] << std::endl;
    


    //"
};
