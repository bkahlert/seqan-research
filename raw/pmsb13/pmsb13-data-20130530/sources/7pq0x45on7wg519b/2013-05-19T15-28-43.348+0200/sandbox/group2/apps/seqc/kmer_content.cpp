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
typedef std::pair<std::string, int> mypair;

struct IntCmp {
    bool operator()(const mypair &lhs, const mypair &rhs) {
        return rhs.second < lhs.second;
    }
};


/**
*Contructor
*Set Constructorvariables minLength and maxLength wich represented the 
*minmal and the maximal length of the kmers
*/
KmerContent::KmerContent(unsigned min, unsigned max)
{
    	minLength = min;
	maxLength = max;
	error=false;
}

/**
*Function for Counting kmers and find their location
*/
int KmerContent::CountKmers(seqan::Dna5String const & seq)
{
	/**
	*check that the kmer length is not higher than the sequenc length
	*/
	if(length(seq)<=minLength){
		error=true;
		return 1;
	}

	unsigned i=minLength;
	/**
	*As long as the kmer length is smaller or equal to the maximal length of the kmer and the sequence length
	*count their observation
	*/	
	while (i<=maxLength&&i<=length(seq)){
		
		/**
		*seqit is an Iterator over the sequence
		*kmerit is an Iterator wich starts at the Position of seqit and goes forward over the seqence until the length of the kmer is reached
		*kmerit helps to save the different kmer strings
		*/
		seqan::Iterator<seqan::Dna5String>::Type seqit=begin(seq);
		seqan::Iterator<seqan::Dna5String>::Type kmerit;
		
		/**
		*As long as seqit don't reachs the end of the sequence 
		*get the kmer at Position of seqit
		*/
        	for(unsigned j=0;j<=length(seq)-i;++j,++seqit){
			kmerit=seqit;
			std::string kmer;

        		for(unsigned kmerrunner=0;kmerrunner<i;++kmerrunner,++kmerit){
				kmer+=*kmerit;
			}
			
			/**
			*Save the whole count of an kmer
			*/
	  		kmerCount[kmer]++;
			
			/**
			*Save the count of an kmer at a specilized position
			*/
	  		kmerPositions[kmer][j]++;
			
		}

	/**
	*count the kmer length one higher
	*/
	++i;
	}

    std::vector<mypair> myvec(kmerCount.begin(), kmerCount.end());
    assert(myvec.size() >= 6);
    std::partial_sort(myvec.begin(), myvec.begin() + 6, myvec.end(), IntCmp());

    for (int i = 0; i < 6; ++i) {
        std::cout << i << ": " << myvec[i].first 
            << "-> " << myvec[i].second << "\n";
    }

return 0;
}


