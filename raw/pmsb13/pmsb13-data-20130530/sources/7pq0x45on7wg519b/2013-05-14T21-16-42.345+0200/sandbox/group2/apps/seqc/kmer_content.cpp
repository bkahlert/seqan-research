#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <cmath>
#include <vector>
#include "kmer_content.h"

KmerContent::KmerContent(unsigned min, unsigned max)
{
    	minLength = min;
	maxLength = max;
}

int KmerContent::CountKmers(seqan::Dna5String const & seq)
{
//seqan::Index<seqan::DnaString, seqan::IndexQGram<seqan::OneGappedShape> > index(seq);
//seqan::stringToShape(seqan::indexShape(index), "1111");

	seqan::Iterator<seqan::Dna5String>::Type it=begin(seq);
	unsigned i=minLength;
	while (i<=maxLength){
        	seqan::Shape<seqan::Dna5> shape;
        	seqan::resize(shape,i);
        	if(length(seq)<=i)
	  	return 1;
        	for(unsigned j=0;j<=length(seq)-i;++j,++it){
	  		seqan::hash(shape,it);
//seqan::hash(seqan::indexShape(index), it);
//for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
//        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;
std::cout << "\n\nShape\n"<< seqan::value(shape) << std::endl;
          //kmerCount[value(shape)]++;
	  //kmerPositions[shape].push_back(j);
//std::cout << "\n\nCount\n"<< kmerCount[shape] << std::endl;
//std::cout << "\n\nPosition\n"<< kmerPositions[shape] << std::endl;
		}
	}
return 0;
}


