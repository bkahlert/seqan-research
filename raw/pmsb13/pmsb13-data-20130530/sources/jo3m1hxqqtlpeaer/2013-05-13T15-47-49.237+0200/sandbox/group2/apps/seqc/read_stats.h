#ifndef __READ_STATS_H
#define __READ_STATS_H

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
struct ReadStats{
    /**
     * 2-dim table that holds the count
     * for each possible score char and 
     * position in read 
     */
    seqan::String<seqan::String <unsigned> > scoreCount;

    /**
     * 2-dim table that holds the count 
     * for each Dna5 element and 
     * position in read
     */
    seqan::String<seqan::String <unsigned> > nucCount;

    /**
     * counter for per read number of gc nuc's
     * position indicates number of gc per read
     * value is the counter itself
     *
     */
    seqan::String<unsigned> gcPerReadCount;

    /**
     * counter for per mean read quality score
     * position indicates score
     * value is the counter itself
     *
     */
    seqan::String<unsigned> perReadQualCount;


    /**
     * holds counts for encountered read lengths, 
     * read lengths as keys 
     */
    std::map <int,int> readLenCount;

    /**
     * The object needs to know expected readlength at construction time,
     * grows when longer sequences are encountered
     */
    
    /**
     * All general job info is kept here. 
     * This map provides job summary data.
     */
    std::map <std::string,std::string> jobParams;
    
    
    /**
     * Contruct an object with an expected/initial readLength,
     * if longer reads are encountered, the object dynamically enlarges 
     */
    ReadStats(unsigned readLength);

    /**
     * collects basic statistics from the reads:
     *  - nucleotide count per position in read for each Dna5 member
     *  - count encountered score chars per position in read
     *  - read length statistics
     *
     */
    void collectReadStats(seqan::CharString const & id, seqan::Dna5String const & seq, seqan::CharString const  & qual);

    /**
     * generates a string representation of collected statistics
     */
    void toString();

    

    void collectQuantiles(seqan::String<unsigned> & cum_intervals, seqan::String<char> & quantiles, unsigned points);

    /**
     * method that accepts reads and queues them for duplication stats
     */
    void collectForDuplicationCheck(seqan::Dna5String & seq);


   /** score offset used */ 
    unsigned scoreOffset;

    /** */
    void grow(unsigned size);

    unsigned getSize();

    private: 
        /** the stats size for registering position-in-read dependent counts, grows with readlength */
        unsigned size;
        
       
        unsigned readCount;
        seqan::StringSet<seqan::Dna5> dupHaystack;
        seqan::Index<seqan::StringSet<seqan::Dna5>, seqan::IndexEsa<> > dupIndex;
        seqan::Finder<seqan::Index<seqan::StringSet<seqan::Dna5>, seqan::IndexEsa<> > > dupFinder;
};

#endif
