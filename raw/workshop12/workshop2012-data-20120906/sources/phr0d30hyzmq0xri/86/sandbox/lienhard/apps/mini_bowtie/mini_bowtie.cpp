// ==========================================================================
//                                mini_bowtie
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.


#include <iostream>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/index_fm.h>
#include <seqan/store.h>

using namespace seqan;
template <typename TIter, typename TStringSet>
int search(TIter & it, TStringSet const & pattern){
  typedef typename Iterator<TStringSet const, Standard>::Type TPatternIter;
  TPatternIter it2;


  return 0;
}

int main(int argc, char *argv[]) 
{
    // type definitions
    typedef String<Dna5> TString;
    typedef StringSet<TString> TStringSet;
    typedef Index<StringSet<TString>, FMIndex<> > TIndex;
    typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIter;

    


    // reading the command line arguments
    if (argc < 3) {
      std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE: mini_bowtie GENOME.fasta READS.fasta OUT.sam" << std::endl;
      return 1;

    }
    StringSet<TString> text;
    
    // declaration and initialization of the fragment store
    FragmentStore<> fragStore;

    // combining the contigs of the reference into one string set
    if (!loadContigs(fragStore, argv[1])) return 1;
    if (!loadReads(fragStore, argv[2])) return 1;

    
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)   
      appendValue(text, fragStore.contigStore[i].seq);
    
    // forward search
    TIndex fmIndex(text);
    TIter it(fmIndex);

    search(it, text.readSeqStore);

    clear(fmIndex);
    clear(it);

    // reversing the sequences for backward search
    reverse(text);
    reverse(fragStore.readSeqStore);


    // backward search
    it = TIter(fmIndex);
    fmIndex = TIndex(text);

    search(it, text.readSeqStore);
    clear(fmIndex);
    clear(it);


    return 0;



}