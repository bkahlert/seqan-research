#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>


ReadStats::ReadStats(unsigned readLength)
{
    seqan::resize(scoreCount, readLength);
    seqan::resize(nucCount,readLength);

        for(unsigned i = 0; i < length(scoreCount); ++i){
            seqan::resize(scoreCount[i], 255, 0); // 2nd dim is ascii chars for scores
            seqan::resize(nucCount[i], 5, 0);     // 2nd dim is elems of Dna5
        }

        std::cout << "Init readStats with readLength " << readLength << std::endl;
};

void ReadStats::toString(){
        std::cout << "== Read Stats == " << std::endl;
        for(unsigned i = 0; i < length(scoreCount); ++i)
            std::cout << scoreCount[i] << std::endl;
};

