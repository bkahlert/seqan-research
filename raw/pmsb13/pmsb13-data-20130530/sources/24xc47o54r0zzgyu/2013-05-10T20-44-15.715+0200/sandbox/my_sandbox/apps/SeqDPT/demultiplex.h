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
/*! \file demultiplex.h
\brief Contains the functions for barcode demultiplexing.
*/

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
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;

typedef Dna5String TAlphabet;

struct seqAtts					//Struct, holding all information to identify a single sequence
{
	unsigned index;
	int barcode;
};

/*!
* Function for controlling the sequences and deleting too short ones. Also checks the barcodes
*/
template <typename TSeqs, typename TIds, typename TBarcodes>
bool check(TSeqs& seqs, TIds& ids, TBarcodes& barcodes) //returns false on errors with the barcodes, true otherwise
{
	unsigned len = length(barcodes[0]);
	for(unsigned i = 1; i < length(barcodes); ++i)
	{
		if(len != length(barcodes[i]))
		{
			std::cerr << "ERROR: Barcodes differ in length. All barcodes must be of equal length.\n";
			return false;	//In this case the programm must be terminated, otherwise erros might occur later or wrong results are likely to be produced
		}
	}
	unsigned limit = length(seqs);
	unsigned i = 0;
	while(i < limit) 
	{
		if(length(seqs[i]) <= len)	
			{
				std::cout << "WARNING: Sequence " << ids[i] << " has length <= barcode length. Sequence will be excluded.\n";
				removeValueById(seqs, i);
				removeValueById(ids, i);
				--limit;
			}
		else ++i;
	}
	return true;
}

template <typename TSeqs, typename TIds, typename Tquals, typename TBarcodes> //Overload for data with quality scores
bool check(TSeqs& seqs, TIds& ids, Tquals& quals, TBarcodes& barcodes) //returns false on errors with the barcodes, true otherwise
{
	unsigned len = length(barcodes[0]);
	for(unsigned i = 1; i < length(barcodes); ++i)
	{
		if(len != length(barcodes[i]))
		{
			std::cerr << "ERROR: Barcodes differ in length. All barcodes must be of equal length.\n";
			return false;	//In this case the programm must be terminated, otherwise erros might occur later or wrong results are likely to be produced
		}
	}
	unsigned limit = length(seqs);
	unsigned i = 0;
	while(i < limit) 
	{
		if(length(seqs[i]) <= len)	
			{
				std::cout << "WARNING: Sequence " << ids[i] << " has length <= barcode length. Sequence will be excluded.\n";
				removeValueById(seqs, i);
				removeValueById(quals, i);
				removeValueById(ids, i);
				--limit;
			}
		else ++i;
	}
	return true;
}

template <typename TSeqs, typename TIds, typename TBarcodes> //Overload for paired-end data
bool check(TSeqs& seqs, TSeqs& seqsRev, TIds& ids, TBarcodes& barcodes) //returns false on errors with the barcodes, true otherwise
{
	unsigned len = length(barcodes[0]);
	for(unsigned i = 1; i < length(barcodes); ++i)
	{
		if(len != length(barcodes[i]))
		{
			std::cerr << "ERROR: Barcodes differ in length. All barcodes must be of equal length.\n";
			return false;	//In this case the programm must be terminated, otherwise erros might occur later or wrong results are likely to be produced
		}
	}
	unsigned limit = length(seqs);
	unsigned i = 0;
	while(i < limit) 
	{
		if(length(seqs[i]) <= len)	
			{
				std::cout << "WARNING: Sequence " << ids[i] << " has length <= barcode length. Sequence will be excluded.\n";
				removeValueById(seqs, i);
				removeValueById(seqsRev, i);
				removeValueById(ids, i);
				--limit;
			}
		else ++i;
	}
	return true;
}

template <typename TSeqs, typename TIds, typename Tquals, typename TBarcodes> //Overload for paired-end data with quality scores
bool check(TSeqs& seqs, TSeqs& seqsRev, TIds& ids, Tquals& quals, TBarcodes& barcodes) //returns false on errors with the barcodes, true otherwise
{
	unsigned len = length(barcodes[0]);
	for(unsigned i = 1; i < length(barcodes); ++i)
	{
		if(len != length(barcodes[i]))
		{
			std::cerr << "ERROR: Barcodes differ in length. All barcodes must be of equal length.\n";
			return false;	//In this case the programm must be terminated, otherwise erros might occur later or wrong results are likely to be produced
		}
	}
	unsigned limit = length(seqs);
	unsigned i = 0;
	while(i < limit) 
	{
		if(length(seqs[i]) <= len)	
			{
				std::cout << "WARNING: Sequence " << ids[i] << " has length <= barcode length. Sequence will be excluded.\n";
				removeValueById(seqs, i);
				removeValueById(seqsRev, i);
				removeValueById(quals, i);
				removeValueById(ids, i);
				--limit;
			}
		else ++i;
	}
	return true;
}

/*!
* Function for extracting the first six bases of all given sequences.
*/
template <typename TSeqs>
TSeqs getPrefix(TSeqs& seqs, unsigned len)
{
	
	TSeqs prefices;
	unsigned limit = length(seqs);
	unsigned i = 0;
	while( i < limit) 
	{
		/*if(len >= length(seqs[i]))	//has to be checked beforehand by seperate function, otherwise might result in faulty mapping of sequences and ids
		{
			std::cout << "WARNING: Sequence " << i << " of current batch has length <= given barcode length. Sequence is beeing excluded." << std::endl;
			removeValueById(seqs, i);
			limit = length(seqs);
		}
		else */
		//{
		appendValue(prefices, prefix(seqs[i], len));
		++i;
		//}
	}
	return prefices;
}

/*!
* Function for exact searching one piece of sequence in the barcode indices
*/
template <typename TSeq, typename TFinder>							//must get the vector/array of finders
int findExactIndex(const TSeq& seq, TFinder& finder)
{
	clear(finder);
	if(find(finder, seq)) 
	{
		return getSeqNo(position(finder));								//returns index of barcode -> ONLY THE FIRST HIT!
	}
	else return -1;														//return -1 if no hit occured
}

/*!
* Function for performing the exact search on all given sequence pieces and barcodes
*/
template <typename TSeqs, typename TFinder>
std::vector<int> findAllExactIndex(const TSeqs& seqs, TFinder& finder) 
{
	std::vector<int> matches(length(seqs));		//initialises the result vector
	for(unsigned i = 0; i < length(seqs); ++i)
	{
		matches[i] = findExactIndex(seqs[i], finder); //Return Vector of barcode indices. Element i belongs to sequence i.
	}
	return matches;
}

/*!
*\brief Function for creating a vector of patterns for approximate searching

The pattern-objects generated use DPSearch and a simple score of 0,-1,-2 (match, mismacht, gap)

\param barcodes a StringSet of barcodes which are DnaStrings
\return A vector of pattern-Objects

*/
template <typename TBarcodes>
std::vector<Pattern<String<Dna>, DPSearch<SimpleScore> > > makePatterns(TBarcodes& barcodes)
{
	std::vector<Pattern<String<Dna>, DPSearch<SimpleScore> > > patterns;
	resize(patterns, length(barcodes));
	for(unsigned i = 0; i < length(barcodes); ++i)
	{
		String<Dna> currentBc = barcodes[i];
		Pattern<String<Dna>, DPSearch<SimpleScore> > pattern(currentBc, SimpleScore(0, -1, -2));
		patterns[i] = pattern;	
	}
	return patterns;
}

/*!
\brief Function for approximate search for one piece of sequence in all barcodes

\param seq Prefix of a sequence.
\param patterns Vector of pattern-objects generated from the barcodes.

\return An integer representing the position of the matched barcode. If no hit occured, -1 is returned.

\remark The given prefix and the barcodes of the pattern-objects must have equal length.

\warning Only the first hit is reported.
*/
template <typename TSeq, typename TPatterns>		//must get the sequence-piece an the vector of patterns(i.e. preprocessed barcodes)
int findApprox(TSeq& seq, TPatterns& patterns)
{
	for (unsigned i = 0; i < length(patterns); ++i)
	{
		Finder<TSeq> finder(seq);			//resets Finder
		if(find(finder, patterns[i], -1))
		{
			return i;						//returns index of matched barcode
		}
	}
	return -1;								//returns -1 if no hit occured
}
		
/*!
\brief Function for approximate search for all pieces of sequence in all barcodes.

\param seqs StringSet of prefices of sequences.
\param patterns Vector of pattern-objects generated from the barcodes.

\return A vector of integers representing the positions of the matched barcodes. 

\remark The n-th element of the vector is associated with the n-th sequence.
\remark The prefices and barcodes of the pattern must have equal length.
\warning Only the first hit of every sequence is reported.

*/
template <typename TSeqs, typename TPatterns>		//must get the StringSet of sequence pieces and the vector of pattern(i.e. preprocessed barcodes)
std::vector<int> findAllApprox(const TSeqs& seqs, TPatterns& patterns )
{
	std::vector<int> matches(length(seqs));			//initialises the result vector
	for(unsigned i = 0; i < length(seqs); ++i)
	{
		matches[i] = findApprox(seqs[i], patterns);
	}
	return matches;
}

/*!
\brief Function for clipping the barcodes from the sequences if they could be matched beforehand

\param seqs StringSet of sequences.
\param matches Result vector of one of the search functions, holding information about the matched barcodes.
\param len unsigned integer representing the length of the barcodes.

\remark The n-th element of seq must be associated with the n-th element of matches.

\warning If the length of the sequences in seq has not been controlled beforehand by check(...) the program might crash.
*/
template <typename TSeqs>			//must get the StringSet of Sequences, the result Vector of a barcode search an the barcode length
void clipBarcodes(TSeqs& seqs, const std::vector<int>& matches, unsigned len)
{
	for(unsigned i = 0; i < length(matches); ++i)
	{
		if(matches[i] != -1)					//only erases barcode from sequence if it could be matched
		{
			erase(seqs[i], 0 , len);
		}
	}	
}

/*!
 \brief overload, if the first len bases of every sequence shall be clipped regardless of the matching(or not matching) barcodes.
 \warning Only to be used if the first len bases of every sequence are sure to hold a barcode.
*/
template <typename TSeqs>			//overload if the first six bases shall be cut in any case
void clipBarcodes(TSeqs& seqs,  int len)
{
	for(unsigned i = 0; i < length(seqs); ++i)
	{
		erase(seqs[i], 0 , len);	
	}	
}

/*!
\brief Function for removing memberless groups from the vector. 

\param groups Two-dimensional vector containing the the seqAtts-objects. Used by the group function.

\remark The vector is not resized. Therefore unused barcodes should be exluded from the input file beforehand to save space an computation time.
*/
void kickEmptyGroups(std::vector<std::vector<seqAtts> >& groups)
{
	unsigned len = length(groups);
	unsigned i = 0;
	while(i < len)
	{
		if(length(groups[i]) == 0)
		{
			groups.erase(groups.begin()+i);		//erases emtpy groups
			--len;								//updates length
		}
		else 
		{
			(groups[i]).shrink_to_fit();		//resizes groups
			++i;
		}
	}
	groups.shrink_to_fit();						//resizes whole vector
}

/*!
\brief Function for sorting the matched and clipped sequences into one vector

The function creates seqAtts-objects from the sequences and barcodes and sorts them into a two dimensional vector.

\param seqs StringSet of the sequences to be grouped
\param ids StringSet of the ids of the sequences
\param matches Vector holding information about the matched barcodes. Produced by one of the searchAll... functions.
\param barcodes StringSet of barcodes.

\return A two-dimensional vector of seqAtts-objects.

\remark The 0-th group (i.e. 0-th collum of the vector) contains the unidentified sequences.
\remark All following groups hold the sequences associated with the i-1-th barcode.

*/
template <typename TSeqs, typename TIds, typename TMatches, typename TBarcodes> //must get the StringSet of sequences (NOT Prefices!), ids, the vector containing the matches and the StringSet of barcodes
std::vector<std::vector<seqAtts> > group(const TSeqs& seqs, const TIds ids, const TMatches& matches, const TBarcodes& barcodes)
{
	std::vector<std::vector<seqAtts> > sortedSequences;			//initialises the result vector
	resize(sortedSequences, length(barcodes)+1);
	for(unsigned i = 0; i < length(matches); ++i)
	{						
		seqAtts attributes;										//initialises struct for attributes of the current sequence...
		attributes.index = i;									//...and assigns the corresponding values
		if(matches[i]!=-1) 
		{
			attributes.barcode = matches[i];					//only assigns barcode if it has been found beforehand...			
		}
		else 
		{
			attributes.barcode = -1;							//...otherwise assigns dummy
		}
		appendValue(sortedSequences[matches[i]+1], attributes);	//adds attributes struct of current sequence to respektive group (offset by 1 is necessary since group 0 is reserved for unidentified sequences)
	}
	kickEmptyGroups(sortedSequences);							//resizes the groups to take minimal space and removes emtpy ones
	return sortedSequences;
}


#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_
