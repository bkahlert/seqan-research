// ==========================================================================
//                                mini_bowtie
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//



#include <iostream>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/index_fm.h>
#include <seqan/store.h>

using namespace seqan;

void search() {};

int main(int argc, char *argv[]) 
{
    // type definitions
    typedef String<Dna5> TString;
    typedef StringSet<TString> TStringSet;
    typedef Index<StringSet<TString>, FMIndex<> > TIndex;
    typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIter;

    // reversing the sequences for backward search
    // backward search
    // reading the command line arguments
    // declaration and initialization of the fragment store
    // forward search
    // combining the contigs of the reference into one string set
    appendValue(text, fragStore.contigStore[i].seq);
    std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE: mini_bowtie GENOME.fasta READS.fasta OUT.sam" << std::endl;
    }
    if (argc < 3) {
    if (!loadContigs(fragStore, argv[1])) return 1;
    if (!loadReads(fragStore, argv[2])) return 1;
    clear(fmIndex);
    clear(fmIndex);
    StringSet<TString> text;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)   
    fmIndex = TIndex(text);
    TIndex fmIndex(text);
    TIter it(fmIndex);
    search();
    search();
    clear(it);
    clear(it);
    reverse(text);
    reverse(fragStore.readSeqStore);
    it = TIter(fmIndex);
    FragmentStore<> fragStore;
    return 0;
    return 1;
}
