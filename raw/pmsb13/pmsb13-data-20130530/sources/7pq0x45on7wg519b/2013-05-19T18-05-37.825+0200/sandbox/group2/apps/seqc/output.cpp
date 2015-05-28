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
#include "kmer_content.h"

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
static const char * SEQC_SEQ_DUPL_DIST_FILE     = "seq_dupl_dist.tsv";
static const char * SEQC_SEQ_DUPL_SEQ_FILE      = "seq_dupl_seq.tsv";
static const char * SEQC_SUMMARY_FILE           = "summary.tsv";
static const char * SEQC_KMER_FILE		= "kmer_content.tsv";
 
// textual conventions
static const char * SEQC_END_MODULE_LN = ">>END_MODULE"; 




/**
 * wrapper: actually writes statistics info to files. 
 */
int TSVWriter::writeAll(ReadStats & stats, KmerContent & kmer, std::string const outdir){
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

    {// GC Content distribution file 
        ofStream.open((outdir + "/" + SEQC_SEQ_GC_CONTENT_FILE).c_str());
        std::stringstream ss;
        writePerReadGcContent(stats, ss);    
        ofStream << ss.str();
        ofStream.close();
    }
    
    {// Overall quality score distribution file 
        ofStream.open((outdir + "/" + SEQC_SEQ_QUAL_DIST_FILE).c_str());
        std::stringstream ss;
        writePerReadQualContent(stats, ss);    
        ofStream << ss.str();
        ofStream.close();
    }
    
    {// Sequence Duplication distribution file 
        ofStream.open((outdir + "/" + SEQC_SEQ_DUPL_DIST_FILE).c_str());
        std::stringstream ss;
        writeSeqDuplDist(stats, ss);    
        ofStream << ss.str();
        ofStream.close();
    }


    {//Sequence Duplication Sequences
        ofStream.open((outdir + "/" + SEQC_SEQ_DUPL_SEQ_FILE).c_str());
        std::stringstream ss;
        writeSeqDuplSeq(stats, ss);    
        ofStream << ss.str();
        ofStream.close();

    }

    {//Kmer Content
        ofStream.open((outdir + "/" + SEQC_KMER_FILE).c_str());
        std::stringstream ss;
        writeKmerContent(kmer, ss);    
        ofStream << ss.str();
        ofStream.close();
    }

    {// Write job summary
        
        ofStream.open((outdir + "/" + SEQC_SUMMARY_FILE).c_str());
        std::stringstream ss;
        writeJobSummary(stats, ss);    
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
    std::map<unsigned, int>::iterator iter; 
  
    for(iter=stats.readLenCount.begin();iter != stats.readLenCount.end(); iter++){
        ss << iter->first << "\t" << iter->second << std::endl;
    }
    

};

void TSVWriter::writePerReadQualContent(ReadStats & stats, std::stringstream & ss){
    ss <<    ">>Per sequence quality scores\tpass";
    ss <<  "#Quality\tCount";


    
};


void TSVWriter::writePerReadGcContent(ReadStats & stats, std::stringstream & ss){
    ss << ">>Per Sequence GC Content\tpass" << std::endl;
    ss << "#GC\tCount" << std::endl;
  
    for(int i=0; i<seqan::length(stats.gcPerReadCount); ++i){
        ss << i << "\t" << stats.gcPerReadCount[i] << std::endl;
    }
};

unsigned sum(seqan::String<unsigned >  counts ){
    unsigned sum = 0; 
    for(unsigned c = 0; c < length(counts);++c){
        sum += counts[c];
    }
    return sum;

};

void TSVWriter::writeNucContent(ReadStats & stats, std::stringstream & gc_content, 
                                std::stringstream & n_content, std::stringstream & nuc_dist )
{
    gc_content << ">>Per base GC content\tpass" << std::endl;
    gc_content << "#Base\t%GC" << std::endl;
    
    n_content << ">>Per base N content\tpass" << std::endl;
    n_content << "#Base\t%GC" << std::endl;

    nuc_dist << ">>Per base sequence content\tpass" << std::endl;
    nuc_dist << "#Base\tG\tA\tT\tC\tN" << std::endl;
        
    // sums for per pos summations
    unsigned sum_a = 0;
    unsigned sum_t = 0;
    unsigned sum_c = 0;
    unsigned sum_g = 0;
    unsigned sum_n = 0;
   
    unsigned sum_all = 0;

    //per position in read:
    for(unsigned pos = 0; pos < length(stats.nucCount); ++pos)
    {

        sum_a = 0;
        sum_t = 0;
        sum_c = 0;
        sum_g = 0;
        sum_n = 0;
        
        sum_c += sum(stats.nucCount[pos][(seqan::Dna5)'C']);
        sum_g += sum(stats.nucCount[pos][(seqan::Dna5)'G']);
        sum_a += sum(stats.nucCount[pos][(seqan::Dna5)'A']);
        sum_t += sum(stats.nucCount[pos][(seqan::Dna5)'T']);
        sum_n += sum(stats.nucCount[pos][(seqan::Dna5)'N']);

        sum_all = sum_c + sum_g + sum_a + sum_t + sum_n;

        if(sum_all > 0)
        {
            gc_content << (pos+1) << "\t"; 
            gc_content << 100*(float)(sum_g+sum_c)/sum_all  << std::endl;

            n_content << (pos+1) << "\t";
            n_content << 100*(float)sum_n/sum_all << std::endl;

            nuc_dist << (pos+1) << "\t";
            nuc_dist 
                << 100*(float)sum_g/sum_all  << "\t"
                << 100*(float)sum_a/sum_all  << "\t"
                << 100*(float)sum_t/sum_all  << "\t"
                << 100*(float)sum_c/sum_all  << "\t"
                << 100*(float)sum_n/sum_all 
                << std::endl;
        }

    }
};



void TSVWriter::writePerPosQual(ReadStats & stats, std::stringstream & ss)
{
    ss << ">>Per base sequence quality\tfail" << std::endl;
    ss << "Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile" 
        << std::endl;
  
            //per position in read:
            for(unsigned pos = 0; pos < length(stats.scoreCount); ++pos){
                // we have to know how many there are counted in total
                // so we can pick specific quantiles later
                __uint64 total_scores = 0;
                
                //we'll cumulate counts here so we have easier access later  
                seqan::String<__uint64> intervals;
               
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
                    ss << (pos+1) << "\t" << ((double)score_sum/total_scores)-stats.scoreOffset ;                
                    ss << "\t" <<  quantiles[9] - stats.scoreOffset;  // median
                    ss << "\t" <<  quantiles[4] - stats.scoreOffset;  // lower quartile
                    ss << "\t" <<  quantiles[14] - stats.scoreOffset;  // upper quartile
                    ss << "\t" <<  quantiles[1] - stats.scoreOffset;  // 10th percentile
                    ss << "\t" <<  quantiles[17] - stats.scoreOffset;  // 90th percentile
                    ss << std::endl;
                }
        }
    
};

void TSVWriter::writeJobSummary(ReadStats & stats, std::stringstream & ss){
        ss << ">>Basic Statistics	pass" << std::endl;
        ss << "#Measure\tValue" << std::endl;
        ss << "Filename\t" << stats.jobParams["Filename"] <<std::endl;	
        ss << "File type\t" << stats.jobParams["File type"] << std::endl;
        
        ss << "Sequence length\t" << stats.readLenCount.begin()->first << std::endl;
        if(stats.readLenCount.begin()->first != stats.readLenCount.rbegin()->first){
            ss << "-" << stats.readLenCount.rbegin()->first; 
        }
        ss << std::endl;
        
        __uint64 sum_at = 0;
        __uint64 sum_gc = 0;
        __uint64 sum_n = 0;

        for(unsigned pos = 0; pos < length(stats.nucCount); ++pos){
            sum_gc += sum(stats.nucCount[pos][(seqan::Dna5)'C']);
            sum_gc += sum(stats.nucCount[pos][(seqan::Dna5)'G']);
            sum_at += sum(stats.nucCount[pos][(seqan::Dna5)'A']);
            sum_at += sum(stats.nucCount[pos][(seqan::Dna5)'T']);
            sum_n += sum(stats.nucCount[pos][(seqan::Dna5)'N']);
        }
        ss << "GC Content\t" << 100 * (double)sum_gc/(sum_gc + sum_at + sum_n) << "%"  << std::endl;
        ss << "N Content\t" << 100 * (double)sum_n/(sum_gc + sum_at + sum_n) << "%"  << std::endl;
        
        ss << "Total Sequences\t" << stats.readCount <<  std::endl;
      
        std::cout << "Calculating Score Range Info ..." << std::endl;
        ScoreInfo score_info;
        guessQualScoreOffset(stats, score_info);
        ss << "Score Range\t " << score_info.description << " [" << (int)score_info.minScore << ":" 
            << (int)score_info.maxScore << "]" << std::endl ;
        

};

/**
 * write Sequence Duplication distribution data
 */
void TSVWriter::writeSeqDuplDist(ReadStats & stats,std::stringstream & ss){
    
    ss << ">>Sequence Duplication Level Distribution\tpass" << std::endl;
    ss << "#Duplication Level\tcount" << std::endl;
    if(stats.duplicates.empty()) 
        return;

    std::map<seqan::Dna5String, unsigned>::iterator iter;
    std::map<unsigned,unsigned> levels;

    for(iter=stats.duplicates.begin();iter != stats.duplicates.end(); iter++){
        ++levels[iter->second];
    }

    for (unsigned l = 1; l<=levels.rbegin()->first;++l){
        ss << l << "\t" << levels[l] << std::endl;
    }
    
}




/**
 * write Sequence Duplication sequences information
 */
void TSVWriter::writeSeqDuplSeq(ReadStats & stats,std::stringstream & ss){
    
    ss << ">>Sequence Duplication Levels\tpass" << std::endl;
    ss << "#Duplicated Sequence\tcount" << std::endl;
    if(stats.duplicates.empty()) 
        return;
    unsigned maxDuplSeqPrint = 50;

    std::map<seqan::Dna5String, unsigned>::iterator iter; 
    typedef std::map<unsigned, seqan::String<seqan::Dna5String> > SMapType;
    SMapType sorted;
    for(iter=stats.duplicates.begin();iter != stats.duplicates.end(); ++iter){
       // ss << iter->first << "\t" << iter->second+1 << std::endl;
        seqan::appendValue(sorted[iter->second],iter->first);
    }
    SMapType::reverse_iterator sit;

    for(sit=sorted.rbegin();sit!=sorted.rend();++sit){
        for(unsigned i = 0; i < length(sit->second)&& (maxDuplSeqPrint > 0); ++i){
            ss << sit->first << "\t" << sit->second[i] << std::endl; 
            --maxDuplSeqPrint;
        }
    }
   
    

}

/**
 * write Kmer Content information
*/
void TSVWriter::writeKmerContent(KmerContent & kmer,std::stringstream & ss){
    
    ss << ">>Kmer Content\tpass" << std::endl;
    ss << "#Position";
    
    std::map<std::string, unsigned>::iterator countIter; 
    std::map<unsigned, std::map<std::string,unsigned> >::iterator kmerPositionIter; 
    std::map<std::string,unsigned>::iterator writeIter; 


    for(countIter=kmer.kmerCount.begin();countIter != kmer.kmerCount.begin(); ++countIter){
        ss<<"\t"<<countIter->first<<std::endl;
    }

    for(kmerPositionIter=kmer.kmerPositions.begin();kmerPositionIter!=kmer.kmerPositions.end();++kmerPositionIter){
	ss<<kmerPositionIter->first;
        for(countIter=kmer.kmerCount.begin();countIter != kmer.kmerCount.begin(); ++countIter){
		
		writeIter=(kmerPositionIter->second).find(countIter->first);
		ss<<"\t"<<writeIter->second;
	}
	ss<<std::endl;

    }
   
}


static const char * _QSCORE_S = "Sanger";
static const char * _QSCORE_L = "Illumnia 1.8+";
static const char * _QSCORE_X = "Solexa";
static const char * _QSCORE_I = "Illumnia 1.3+";


/**
 * use existing quality scores to guess the scoring offset
 * result 
 */
void TSVWriter::guessQualScoreOffset(ReadStats & stats, ScoreInfo & score_info){
    score_info.minScore = 255;
    score_info.maxScore = 0;
   
    //Scoring system intervals
    char intervals[4][3] = {
        {'S',33,73},
        {'L',33,74},
        {'X',59,104},
        {'I',64,104}
    };

    std::map<char,seqan::String<char> > names;
    names['S'] = "Sanger";
    names['L'] = "Illumina 1.8+";
    names['X'] = "Solexa";
    names['I'] = "Illumina 1.3+";
    // find min and max scores
    for(unsigned s = 0; s < 255; ++s){
        for(unsigned pos = 0; pos < stats.getSize(); ++pos){
            if(stats.scoreCount[pos][s]>0){
                if(s < score_info.minScore){
                    score_info.minScore = s;
                }
                if(s > score_info.maxScore){
                    score_info.maxScore = s;
                }
            }
        }
    }

    // determining offset
    if(score_info.minScore < 59){
        score_info.offset = 33;
    } else {
        score_info.offset = 64; 
    }

    std::string possible; 
    for(int i=0; i<(sizeof(intervals)/sizeof(*intervals));++i){
        if(intervals[i][1] <= score_info.minScore && intervals[i][2] >= score_info.maxScore){
            seqan::append(possible,", ");
            seqan::append(possible,names[(char)intervals[i][0]]);
        }
    }
    score_info.description = possible.substr(2,std::string::npos);

    


};

