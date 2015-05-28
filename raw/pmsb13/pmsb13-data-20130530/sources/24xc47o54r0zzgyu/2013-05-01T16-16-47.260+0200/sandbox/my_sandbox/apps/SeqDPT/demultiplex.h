// ==========================================================================
//                               readTrimming.h
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
// Author: Sebastian Roskosch 
// ==========================================================================

#ifndef SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;

typedef Dna5String TAlphabet;


//template <typename TBcs, typename TIndices, typename TFinders> 
//void makeExactIndices(const TBcs& barcodes, TIndices& indices, TFinders& finders)
//{
//	Index<TAlphabet,IndexEsa<> > indices(barcodes);
//	Finder<esaIndex::Type, IndexEsa<> > finders(indices);
//}

/*!
* Function for searching 1 piece of sequence in the barcode indices
*/
template <typename TSeq, typename TFinders>							//must get the vector/array of finders
int findExactIndex(const TSeq& seq, TFinders& finders)
{
	clear(finders);
	if(find(finders, seq)) 
	{
		return getSeqNo(position(finders));								//returns index of barcode -> ONLY THE FIRST HIT!
	}
	else return -1;														//return -1 if no hit occured
}

/*!
* Function for performing the search on all given sequence pieces and indices
*/
template <typename TSeqs, typename TFinder>
std::vector<std::pair<int, int> > findAllExactIndex(const TSeqs& seqs, TFinder& finder) 
{
	std::vector<std::pair<int, int> > matches(length(seqs));		//initialises the result vector
	//Iterator<TSeqs, Rooted>::Type it = begin(seqs);
	/*for(goBegin(it); !atEnd(it); goNext(it))
	{
		appendValue(matches, std::make_pair(position(it), findExactIndex(it, finder)));
	}*/
	return matches;
}


#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_
