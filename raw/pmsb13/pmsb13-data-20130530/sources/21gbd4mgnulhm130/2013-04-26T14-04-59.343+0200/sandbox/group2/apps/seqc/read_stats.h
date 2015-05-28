#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

struct ReadStats{
    
    seqan::String<seqan::String <unsigned> > scoreCount;
    seqan::String<seqan::String <unsigned> > nucCount;

    ReadStats(unsigned readLength);
    void collectReadStats(seqan::CharString const & id, seqan::Dna5String const & seq, seqan::CharString const  & qual);

    void toString();
};

