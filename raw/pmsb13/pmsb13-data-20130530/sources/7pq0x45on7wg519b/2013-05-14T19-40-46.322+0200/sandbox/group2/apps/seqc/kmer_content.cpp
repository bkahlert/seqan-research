#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <cmath>
#include "kmer_content.h"

KmerContent::KmerContent(unsigned min, unsigned max)
{
    	minLength = min;
	maxLength = max;
};

void KmerContent::CountKmers(seqan::Dna5String const & seq)
{
    seqan::Iterator<DnaString>::Type it=begin(seq);
    unsigned i=min;
    while (i<=max){
        seqan::Shape<seqan::Dna5> shape;
        seqan::resize(shape,i);
        if(length(seq)<=i)
	  return 0;
        for(unsigned j=0;j<=length(seq)-i;++j,++it){
	  seqan::hash(shape,it);
std::cout << "\n\nShape\n"<< shape << std::endl;
          //kmerCount[shape]++;
	  kmerPositions[shape].push_back(j);
std::cout << "\n\nCount\n"<< kmerCount[shape] << std::endl;
std::cout << "\n\nPosition\n"<< kmerPositions[shape] << std::endl;
	}
    }

};


