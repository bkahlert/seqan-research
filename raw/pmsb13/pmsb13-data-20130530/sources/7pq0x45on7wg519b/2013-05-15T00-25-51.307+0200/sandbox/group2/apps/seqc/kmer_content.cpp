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
	if(length(seq)<=minLength){
		return 1;
}

	seqan::Index<seqan::DnaString, seqan::IndexQGram<seqan::OneGappedShape> > index(seq);
	
	std::string bitmap;
//	std::string kmer;

	unsigned i=minLength;
	while (i<=maxLength&&i<=length(seq)){
		for(unsigned bitmapper=0;bitmapper<i;++bitmapper){
			bitmap+="1";
		}
		seqan::stringToShape(seqan::indexShape(index), bitmap);
//std::cout << "bitmap\t"<< bitmap <<"\n\n"<< std::endl;
		seqan::Iterator<seqan::Dna5String>::Type seqit=begin(seq);
		seqan::Iterator<seqan::Dna5String>::Type kmerit;
//std::cout <<seq<< std::endl;
        	seqan::Shape<seqan::Dna5> shape;
        	seqan::resize(shape,i);

        	for(unsigned j=0;j<=length(seq)-(i+1);++j,++seqit){
		kmerit=seqit;
		std::string kmer;
        		for(unsigned kmerrunner=0;kmerrunner<i;++kmerrunner,++kmerit){
				kmer+=*kmerit;
			}
//std::cout << "\n\nbefor it_val\t"<< kmer << std::endl;
	  		seqan::hash(shape,seqit);
//std::cout << "after it_val\t"<< *it << std::endl;
//seqan::shapeToString(bitmap,shape);
//seqan::hash(seqan::indexShape(index), "ACC");

//for (unsigned x = 0; x < length(getOccurrences(index, indexShape(index))); ++x){
//std::cout << "\n\nbla\n"<< std::endl;
       //std::cout << getOccurrences(index, indexShape(index))[x] << std::endl;}
//std::cout << "\n\nShape\n"<< seqan::value(shape)<< std::endl;
          //kmerCount[value(shape)]++;
	  		kmerPositions[kmer].push_back(j);
		}
//std::cout << "\n\nCount\n"<< kmerCount[shape] << std::endl;
//std::cout << "\n\nPosition\n"<< kmerPositions[shape] << std::endl;
for(std::map<std::string,std::vector <int> >::const_iterator it = kmerPositions.begin(); it != kmerPositions.end(); ++it) 
{
    std::cout << "coord: " << it->first <<std::endl;
//for(unsigned ixi = 0; ixi<(it->second).length();++ixi)
  // std::cout << "property: " << it->second[ixi];
//}
vector<int> vInnererVektor = *it->second;
vector<int>::iterator itInnen;    
    for( itInnen = vInnererVektor.begin(); itInnen != vInnererVektor.end(); itInnen++)
    {
        int iValue = *itInnen;
        std::cout << "property: " << iValue;
        // tue irgendwas mit iValue
    }
++i;
	}
return 0;
}


