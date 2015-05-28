// ==========================================================================
//                              bowtie_tutorial
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================

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
    typedef Index<StringSet<TString>, FmIndex<> > TIndex;
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
	
	TIndex fmIndex(text);
	TIter it(fmIndex);

        // forward search
	search();
    	clear(it);
	clear(fmIndex);
	// reversing the sequences for backward search
	reverse(fragStore.readSeqStore);
	reverse(text);

	fmIndex = TIndex(text);
	it = TIter(fmIndex);
	// backward search
	search();
        clear(it);
	clear(fmIndex);

    return 0;
   
}