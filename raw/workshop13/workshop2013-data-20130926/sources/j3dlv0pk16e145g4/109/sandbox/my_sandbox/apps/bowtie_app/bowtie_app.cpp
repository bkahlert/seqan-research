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
void search(TIter & it, TStringSet const & pattern)
{
    for (unsigned i=0; i < length(pattern); ++i)
    {
        unsigned j=length(pattern[i]-1);
        while(goDown(it, pattern[i][j]) && j>= length(pattern)/2)
        {
            --j;
        }
        goRoot(it);
    }
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
    // declaration and initialization of the fragment store
    FragmentStore<> fragStore;
    if (!loadContigs(fragStore, argv[1])) return 1;
    if (!loadReads(fragStore, argv[2])) return 1;
    // combining the contigs of the reference into one string set
    StringSet<TString> text;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)   
        appendValue(text, fragStore.contigStore[i].seq);
    // forward search
    TIndex fmIndex(text);
    TIter it(fmIndex);
    search();
    // reversing the sequences for backward search
    reverse(text);
    reverse(fragStore.readSeqStore);
    // backward search
    clear(fmIndex);
    fmIndex = TIndex(text);
    clear(it);
    it = TIter(fmIndex);
    search();
    clear(fmIndex);
    fmIndex = TIndex(text);
    clear(it);
    it = TIter(fmIndex);
    clear(fmIndex);
    clear(it);

    return 0;

}