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
            seqan::resize(nucCount[i], 5, 0);     // 2nd dim is elems of Dna5
        }

  //      std::cout << "Init readStats with readLength " << readLength << std::endl;
};

void ReadStats::collectReadStats(seqan::CharString const & id, seqan::Dna5String const & seq, seqan::CharString const  & qual)
{
    SEQAN_ASSERT_EQ(length(seq), length(qual), "read sequence length and quality String length must match");
    for(unsigned pos = 0; pos < length(seq); ++pos){
        ++scoreCount[pos][qual[pos]];
        ++nucCount[pos][seq[pos]];
        std::cout << "In pos ["<< pos << "] " <<  seq[pos] << ":" << qual[pos] << std::endl;
    }
}

void ReadStats::toString(){
        std::cout << "== Read Stats == " << std::endl;
        for(unsigned i = 0; i < length(scoreCount); ++i)
            std::cout << scoreCount[i] << std::endl;
};

