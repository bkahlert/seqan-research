#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <cmath>
#include "read_stats.h"

ReadStats::ReadStats(unsigned readLength)
{
    size = 0;
    readCount = 0;

    //grow datastructures using positional counters 
    grow(readLength);

    scoreOffset = 33;

    // for storing the mean read quality counts
    seqan::resize(perReadQualCount, 255,0);
    
    //initialize empty index
    seqan::indexText(dupIndex);

    //initialize default duplicates haystack size
    dupHaystackSize = 200000;
    
}

void ReadStats::collectReadStats(seqan::CharString const & id, seqan::Dna5String const & seq, seqan::CharString const  & qual)
{
    SEQAN_ASSERT_EQ(length(seq), length(qual));
    id;
    unsigned thelen = length(qual);
    // allow for dynamic growth of the internal structures when we encounter 
    // sequences longer than expected
    if(size<thelen){
        grow(thelen);
    }
    SEQAN_ASSERT_GEQ(size, length(qual));

    // update readLength counter
    readLenCount[thelen]++;

    // update positional counters 
    unsigned gc_per_read = 0;
    unsigned score_sum = 0;
    for(unsigned pos = 0; pos < thelen; ++pos){
        ++scoreCount[pos][qual[pos]];
        ++nucCount[pos][seqan::ordValue(seq[pos])];
        if((seq[pos] == (seqan::Dna5)'G') || (seq[pos] == (seqan::Dna5)'C')){
            ++gc_per_read;
        }
        score_sum += qual[pos];
    }

    //update per read counter for gc content
    ++gcPerReadCount[gc_per_read];
    
    //update per read counter for mean qual score
    //round mean to closest int
    unsigned mean =  std::floor(0.5 + ((double)score_sum/thelen));
    ++perReadQualCount[mean];
}




void ReadStats::toString(){
        std::cout << "== Read Stats == " << std::endl;
        for(unsigned i = 0; i < length(scoreCount); ++i)
            std::cout << scoreCount[i] << std::endl;
}

void ReadStats::collectQuantiles(seqan::String<__uint64> & cum_intervals, seqan::String<char> & quantiles, unsigned points)
{
    seqan::resize(quantiles, points);
    
    //we assume the sum of all counts is the last elem in cum_intervals
    __uint64 total_sum = cum_intervals[length(cum_intervals)-1];
   
    //steps to walk: points, the last interval can be as large as 2(total_sum / points )-1
    seqan::String<unsigned> quantities;
    seqan::resize(quantities, points);

    for(unsigned i = 0; i<points; ++i){
        quantities[i] = ((float)(i+1)/points)*total_sum;
    }
    
    //search for each quantity in cum_intervals: 
    //choose the first bucket in cum_intervals that has a larger count than 
    unsigned quant_find = 0;
    char s = 63;
    while( (unsigned)s<length(cum_intervals) && quant_find < length(quantities)){
        if(quantities[quant_find] <= cum_intervals[s]){
            quantiles[quant_find] = s;
            ++quant_find;
        } else {
            ++s;
        }
    }

}


void ReadStats::grow(unsigned newsize){
    if (newsize <= size){
        return;
    }
    
    seqan::resize(scoreCount, newsize);
    seqan::resize(nucCount,newsize);
       
    for(unsigned i = size; i < newsize; ++i){
            // 2nd dim is ascii chars for scores
            seqan::resize(scoreCount[i], 255, 0);             
           
            // 2nd dim is elems of Dna5
            seqan::resize(nucCount[i], seqan::ValueSize<seqan::Dna5>::VALUE, 0);
        }

    // for gc count per read, initialize with 0 count, 
    // can anyting in[0,readlength]
    // so we need readlength+1 counters
    seqan::resize(gcPerReadCount,newsize+1, 0);
   

    size = newsize;
}

unsigned ReadStats::getSize()
{
    return size;
}

/**
 * find duplicates
 * accepts reads to check for duplicates and accumulates 
 * duplicate counts in ReadStats::duplicates
 * uses the first ReadStats::dupHaystacksize reads as
 * haystack
 */
void ReadStats::dupCheck(seqan::Dna5String const & seq){
    
    if(readCount < dupHaystackSize){
        seqan::appendValue(dupHaystack, seq);
        seqan::appendValue(seqan::indexText(dupIndex), seq);
    } else if(readCount == dupHaystackSize){
        seqan::appendValue(seqan::indexText(dupIndex), seq);
        //create the index
        seqan::indexCreate(dupIndex,seqan::EsaSA());
        //check for duplicates in the haystack itself
        //must find exactly once
        
        seqan::Finder<seqan::Index<seqan::StringSet<seqan::Dna5String>, seqan::IndexEsa<> > > dupFinder(dupIndex);
        for(unsigned i=0;i<seqan::length(seqan::indexText(dupIndex)); ++i){
            seqan::clear(dupFinder);
            
            // find first occurence that we know about
            seqan::find(dupFinder,seqan::indexText(dupIndex)[i]);

            // if there are others, add them to the duplication list
            if(seqan::find(dupFinder,seqan::indexText(dupIndex)[i])){
                ++duplicates[seqan::indexText(dupIndex)[i]];
            }
        }
    } else {
    
        seqan::Finder<seqan::Index<seqan::StringSet<seqan::Dna5String>, seqan::IndexEsa<> > > dupFinder(dupIndex);
    
        if(seqan::find(dupFinder, seq)){
            ++duplicates[seq];    
        }
    }
    
    ++readCount;

    
}


