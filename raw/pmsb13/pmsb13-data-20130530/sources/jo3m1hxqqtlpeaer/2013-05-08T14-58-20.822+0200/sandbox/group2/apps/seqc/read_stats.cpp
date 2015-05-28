#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <cmath>
#include "read_stats.h"

ReadStats::ReadStats(unsigned readLength)
{
    seqan::resize(scoreCount, readLength);
    seqan::resize(nucCount,readLength);
    scoreOffset = 33;
        for(unsigned i = 0; i < length(scoreCount); ++i){
            // 2nd dim is ascii chars for scores
            seqan::resize(scoreCount[i], 255, 0);             
           
            // 2nd dim is elems of Dna5
            seqan::resize(nucCount[i], seqan::ValueSize<seqan::Dna5>::VALUE, 0);
        }

    // for gc count per read, initialize with 0 count, 
    // can anyting in[0,readlength]
    // so we need readlength+1 counters
        seqan::resize(gcPerReadCount,readLength+1, 0);
        
        // for storing the mean read quality counts
        seqan::resize(perReadQualCount, 255,0);
};

void ReadStats::collectReadStats(seqan::CharString const & id, seqan::Dna5String const & seq, seqan::CharString const  & qual)
{
    SEQAN_ASSERT_EQ(length(seq), length(qual));
    readLenCount[length(seq)]++;
    unsigned gc_per_read = 0;
    unsigned score_sum = 0;
    for(unsigned pos = 0; pos < length(seq); ++pos){
        ++scoreCount[pos][qual[pos]];
        ++nucCount[pos][seqan::ordValue(seq[pos])];
        if((seq[pos] == (seqan::Dna5)'G') || (seq[pos] == (seqan::Dna5)'C')){
            ++gc_per_read;
        }
        score_sum += qual[pos];
    }
    ++gcPerReadCount[gc_per_read];
    //round mean to closest int
    unsigned mean = std::floor(0.5 + ((double)score_sum/length(qual)));
    ++perReadQualCount[mean];
};




void ReadStats::toString(){
        std::cout << "== Read Stats == " << std::endl;
        for(unsigned i = 0; i < length(scoreCount); ++i)
            std::cout << scoreCount[i] << std::endl;
};

void ReadStats::collectQuantiles(seqan::String<unsigned> & cum_intervals, seqan::String<char> & quantiles, unsigned points)
{
    seqan::resize(quantiles, points);
    
    //we assume the sum of all counts is the last elem in cum_intervals
    unsigned total_sum = cum_intervals[length(cum_intervals)-1];
   
    //steps to walk: points, the last interval can be as large as 2(total_sum / points )-1
    seqan::String<unsigned> quantities;
    seqan::resize(quantities, points);

    for(unsigned i = 0; i<points; ++i){
        quantities[i] = ((float)(i+1)/points)*total_sum;
    }
    
    //search for each quantity in cum_intervals: 
    //choose the first bucket in cum_intervals that has a larger count than 
    unsigned quant_find = 0;
    char s = 63;
    while( s<length(cum_intervals) && quant_find < length(quantities)){
        if(quantities[quant_find] <= cum_intervals[s]){
            quantiles[quant_find] = s;
            ++quant_find;
        } else {
            ++s;
        }
    }

};


void ReadStats::grow(unsigned size){


};
