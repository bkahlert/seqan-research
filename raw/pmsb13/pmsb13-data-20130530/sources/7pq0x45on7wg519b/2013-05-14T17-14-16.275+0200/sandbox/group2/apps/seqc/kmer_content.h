#ifndef __KMER_CONTENT_H
#define __KMER_CONTENT_H

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <map>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
/** 
 * Collects and holds all information gathered from input
 *
 *
 *
 */
struct KmerContent{
    	
    /**
     * holds counts for encountered kmers, 
     * shape as keys 
     */
    std::map <seqan::Shape<seqan::Dna5>,int> kmerCount;

    /**
     * The object needs to know the minimal and maximal kmer-length at construction time
     */
    KmerContent(unsigned min, unsigned max);

    void CountKmers(seqan::Dna5String const & seq);


    std::map <seqan::Shape<seqan::Dna5>,std::vector<int> > kmerPositions;



    private: 
        /** the minimal kmer length for searching the kmers in all reads*/
        unsigned minLength;
	/** the maximal kmer length for searching the kmers in all reads*/
	unsigned maxLength;
	
};

#endif
