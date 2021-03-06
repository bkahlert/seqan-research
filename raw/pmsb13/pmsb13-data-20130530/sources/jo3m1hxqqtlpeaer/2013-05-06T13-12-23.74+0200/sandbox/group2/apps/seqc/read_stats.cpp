#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "read_stats.h"

ReadStats::ReadStats(unsigned readLength)
{
    seqan::resize(scoreCount, readLength);
    seqan::resize(nucCount,readLength);

        for(unsigned i = 0; i < length(scoreCount); ++i){
            seqan::resize(scoreCount[i], 255, 0); // 2nd dim is ascii chars for scores
            seqan::resize(nucCount[i], seqan::ValueSize<seqan::Dna5>::VALUE, 0);     // 2nd dim is elems of Dna5
        }

  //      std::cout << "Init readStats with readLength " << readLength << std::endl;
};

void ReadStats::collectReadStats(seqan::CharString const & id, seqan::Dna5String const & seq, seqan::CharString const  & qual)
{
    SEQAN_ASSERT_EQ(length(seq), length(qual));
    readLenCount[length(seq)]++;
    for(unsigned pos = 0; pos < length(seq); ++pos){
        ++scoreCount[pos][qual[pos]];
        ++nucCount[pos][seqan::ordValue(seq[pos])];
      //  std::cout << "In pos ["<< pos << "] " <<  seq[pos] << ":" << qual[pos] << std::endl;
    }
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
   
    std::cout << "Working with total sum: " << total_sum << std::endl;

    //steps to walk: points, the last interval can be as large as 2(total_sum / points )-1
    seqan::String<unsigned> quantities;
    seqan::resize(quantities, points+1);

    std::cout << "Searching for the following points in distribution: " << std::endl;
    for(unsigned i = 0; i<points; ++i){
        quantities[i] = ((float)(i+1)/points)*total_sum;
        std::cout << quantities[i] << " " << ((float)(i+1)/points)*total_sum << std::endl;
    }

    //search for each quantity in cum_intervals: 
    //choose the first bucket in cum_intervals that has a larger count than 
    unsigned quant_find = 0;
    for(unsigned s = 0; s<length(cum_intervals) && quant_find<length(quantities) ;++s){
        std::cout << "searching for " << quantities[quant_find] << " in quantile no " << s << ": " << quantiles[s] << std::endl;
        if(quantities[quant_find] <= quantiles[s]){
            std::cout << " found quantities [" << quant_find << ":" << quantities[quant_find] 
                << "] in cum intervals bucket [" << s << ": " << cum_intervals[s] << "]" << std::endl;
            ++quant_find;
            --s;

        }
    }

};
