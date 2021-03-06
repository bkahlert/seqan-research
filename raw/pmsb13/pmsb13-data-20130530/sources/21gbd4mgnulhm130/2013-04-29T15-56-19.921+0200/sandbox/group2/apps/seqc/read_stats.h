#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <map>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
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
     * holds counts for encountered read lengths, 
     * read lengths as keys 
     */
    std::map <int,int> readLenCount;

    /**
     * The object needs to know maximal readlength at construction time 
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
};

