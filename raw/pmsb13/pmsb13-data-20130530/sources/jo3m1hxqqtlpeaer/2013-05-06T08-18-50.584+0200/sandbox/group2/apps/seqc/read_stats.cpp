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
    seqan::resize(quantities, points);

    for(unsigned i = 0; i<points; ++i){
        quantities[i] = (i/points)*total_sum;
        std::cout << quantities[i] << " " << (i/points)*total_sum << std::endl;
    }

};
