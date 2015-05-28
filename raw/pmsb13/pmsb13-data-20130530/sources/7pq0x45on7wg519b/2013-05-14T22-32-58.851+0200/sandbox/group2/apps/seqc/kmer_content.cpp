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
seqan::Index<seqan::DnaString, seqan::IndexQGram<seqan::OneGappedShape> > index(seq);
seqan::stringToShape(seqan::indexShape(index), "1111");
//std::string bitmap;

	unsigned i=minLength;
	while (i<=maxLength){
	seqan::Iterator<seqan::Dna5String>::Type it=begin(seq);
std::cout << "\n\nhallo\n"<< std::endl;
        	seqan::Shape<seqan::Dna5> shape;
        	seqan::resize(shape,i);
        	if(length(seq)<=i)
	  	return 1;
        	for(unsigned j=0;j<=length(seq)-(i+1);++j,++it){
//	  		seqan::hash(shape,it);
//seqan::shapeToString(bitmap,shape);
seqan::hash(seqan::indexShape(index), &it);
std::cout << "\n\nit_val\n"<< *it << std::endl;
for (unsigned x = 0; x < length(getOccurrences(index, indexShape(index))); ++x){
std::cout << "\n\nbla\n"<< std::endl;
       std::cout << getOccurrences(index, indexShape(index))[x] << std::endl;}
//std::cout << "\n\nShape\n"<< bitmap << std::endl;
          //kmerCount[value(shape)]++;
	  //kmerPositions[shape].push_back(j);
//std::cout << "\n\nCount\n"<< kmerCount[shape] << std::endl;
//std::cout << "\n\nPosition\n"<< kmerPositions[shape] << std::endl;
		}
++i;
	}
return 0;
}


