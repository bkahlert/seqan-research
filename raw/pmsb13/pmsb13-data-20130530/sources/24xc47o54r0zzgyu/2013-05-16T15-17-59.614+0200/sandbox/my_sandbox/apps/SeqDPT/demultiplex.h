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
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;

typedef String<Dna5Q> TAlphabet;

// ============================================================================
// Functions
// ============================================================================

/*!
\brief Function for controlling the sequences and deleting too short ones. Also checks the barcodes

The Function deletes the ID and sequence from the given StringSet if the sequence length is <= barcode length.

\param seqs StringSet of forward reads.
\param seqsRev StringSet of backward reads.
\param ids StringSet of IDs.
\param barcodes StringSet of barcodes
\return A bool: \b false on different lengths of the barcodes, \b true otherwise.
\remark The n-th entries of all given StringSets must be associated.
\warning If the programm is not terminated uppon return of \b false, crashes or wrong results are to be exspected.  
*/
template <typename TSeqs, typename TIds, typename TBarcodes> 
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

/*!
\brief Overload of ::check(TSeqs& seqs, TSeqs& seqsRev, TIds& ids, TBarcodes& barcodes) for single-end data.
\param seqs StringSet of reads.
\param ids StringSet of IDs.
\param barcodes StringSet of barcodes
*/
template <typename TSeqs, typename TIds, typename TBarcodes> 
bool check(TSeqs& seqs, TIds& ids,TBarcodes& barcodes) 
{
	unsigned len = length(barcodes[0]);
	for(unsigned i = 1; i < length(barcodes); ++i)
	{
		if(len != length(barcodes[i]))
		{
			std::cerr << "ERROR: Barcodes differ in length. All barcodes must be of equal length.\n";
			return false;
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


/*!
\brief Function for extracting the prefices of all given sequences.
\param seqs StringSet of sequences.
\param len Unsigned integer representing the desired prefix length.
\return A StringSet of prefices of length \b len
*/
template <typename TSeqs>
TSeqs getPrefix(TSeqs& seqs, unsigned len)
{
	TSeqs prefices;
	unsigned limit = length(seqs);
	unsigned i = 0;
	while( i < limit) 
	{
		appendValue(prefices, prefix(seqs[i], len));
		++i;
	}
	return prefices;
}

/*!
\brief Function for exact search for one prefix in the barcode indices.
\param prefix Prefix of a sequence.
\param finder Finder-object holding the preprocessed barcodes.
\return An integer representing the position of the matching barcode. If no hit occured \b -1 is returned.
\warning Only the first hit is reported.
*/
template <typename TPrefix, typename TFinder>							
int findExactIndex(const TPrefix& prefix, TFinder& finder)
{
	clear(finder);
	if(find(finder, prefix)) 
	{
		return getSeqNo(position(finder));								//returns index of barcode -> ONLY THE FIRST HIT!
	}
	else return -1;														//return -1 if no hit occured
}

/*!
\brief Function for performing ::findAllExactIndex on all given prefices and barcodes
\param prefices StringSet of prefices.
\param finder Finder-object holding the preprocessed barcodes.
\return A Vector of integers representing the positions of the matched barcodes. 
If a sequence could not be matched, the position is replaced by \b -1.
\remark The n-th element of the resulting vector is associated with the n-th prefix/sequence.
\warning Only the first hit of each prefix is repoted.
*/
template <typename TPrefices, typename TFinder>
std::vector<int> findAllExactIndex(const TPrefices& prefices, TFinder& finder) 
{
	std::vector<int> matches(length(prefices));				//initialises the result vector
	for(unsigned i = 0; i < length(prefices); ++i)
	{
		matches[i] = findExactIndex(prefices[i], finder);	//returns vector of barcode indices. Element i belongs to sequence i.
	}
	return matches;
}

/*!
*\brief Function for creating a vector of patterns used by ::findApprox and ::findApproxAll.

The pattern-objects generated use DPSearch and a simple score of 0,-1,-2 (match, mismatch, gap)
\param barcodes StringSet of barcodes.
\return A vector of pattern-Objects.
*/
template <typename TBarcodes>
std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > makePatterns(TBarcodes& barcodes)
{
	std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > patterns;
	resize(patterns, length(barcodes));
	for(unsigned i = 0; i < length(barcodes); ++i)
	{
		String<Dna5Q> currentBc = barcodes[i];
		Pattern<String<Dna5Q>, DPSearch<SimpleScore> > pattern(currentBc, SimpleScore(0, -1, -2));
		patterns[i] = pattern;	
	}
	return patterns;
}

/*!
\brief Function for approximate search for one prefix in all barcodes.
\param prefix Prefix of a sequence.
\param patterns Vector of pattern-objects generated by ::makePatterns.
\return An integer representing the position of the matched barcode. If no hit occured \b -1 is returned.
\warning Only the first hit is reported.
*/
template <typename TPrefix, typename TPatterns>
int findApprox(TPrefix& prefix, TPatterns& patterns)
{
	for (unsigned i = 0; i < length(patterns); ++i)
	{
		Finder<TPrefix> finder(prefix);			//resets Finder
		if(find(finder, patterns[i], -1))
		{
			return i;						//returns index of matched barcode
		}
	}
	return -1;								//returns -1 if no hit occured
}
		
/*!
\brief Function for performing ::findApprox on all given prefices and barcodes.
\param prefices StringSet of prefices.
\param patterns Vector of pattern-objects generated by ::makePatterns.
\return A vector of integers representing the positions of the matched barcodes. 
\remark The n-th element of the vector is associated with the n-th sequence.
\warning Only the first hit of each is reported.
*/
template <typename TPrefices, typename TPatterns>		//must get the StringSet of prefices and the vector of patterns(i.e. preprocessed barcodes)
std::vector<int> findAllApprox(const TPrefices& prefices, TPatterns& patterns )
{
	std::vector<int> matches(length(prefices));			//initialises the result vector
	for(unsigned i = 0; i < length(prefices); ++i)
	{
		matches[i] = findApprox(prefices[i], patterns);
	}
	return matches;
}

/*!
\brief Function for clipping the barcodes from the sequences if they could be matched beforehand.
\param seqs StringSet of sequences.
\param matches Vector produced by ::findAllExactIndex or ::findAllApprox, holding information about the matched barcodes.
\param len Unsigned integer representing the length of the barcodes.
\remark The n-th element of seq must be associated with the n-th element of matches.
\warning If the length of the sequences in seq has not been controlled beforehand by ::check the program might crash.
*/
template <typename TSeqs>			
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
\brief Overload of ::clipBarcodes<seqs, matches, len> if the first len bases of every sequence shall be clipped regardless of the matching(or not matching) barcodes.
\param seqs StringSet of sequences.
\param len unsigned integer representing the length of the barcodes.
\warning Only to be used if the first len bases of every sequence are sure to hold a barcode.
\warning If the length of the sequences in seq has not been controlled beforehand by ::check the program might crash.
*/
template <typename TSeqs>			
void clipBarcodes(TSeqs& seqs,  int len)
{
	for(unsigned i = 0; i < length(seqs); ++i)
	{
		erase(seqs[i], 0 , len);	
	}	
}

/*!
\brief Function for resizing the groups of the vector produced by ::group and used by it. 
\param groups Two-dimensional vector containing integers representing the indices of the sequences. Produced by ::group.
*/
void resizeGroups(std::vector<std::vector<int> >& groups)
{
	for(unsigned i = 0; i < length(groups); ++i)
	{
		(groups[i]).shrink_to_fit();			//resizes groups
	}
	groups.shrink_to_fit();						//resizes the whole vector
}

/*!
\brief Function for sorting the matched and clipped sequences into one vector.

The function takes the results of ::findAllExactIndex or ::findAllApprox and the sequnce indices into a vector

\param matches Vector holding information about the matched barcodes. Produced by ::findAllExactIndex or ::findAllApprox functions.
\param barcodes StringSet of barcodes.
\return A two-dimensional vector of integers representing the indices of the sequences.
\remark The 0-th group (i.e. 0-th collum of the vector) contains the unidentified sequences.
\remark All following groups hold the sequences associated with the i-1-th barcode.
*/
template <typename TSeqs, typename TMatches, typename TBarcodes> //must get the StringSet of sequences (NOT Prefices!), the vector containing the matches and the StringSet of barcodes
std::vector<std::vector<int> > group(TSeqs& seqs, const TMatches& matches, const TBarcodes& barcodes)
{
	std::vector<std::vector<int> > sortedSequences;			//initialises the result vector
	resize(sortedSequences, length(barcodes)+1);
	for(unsigned i = 0; i < length(matches); ++i)
	{						
		appendValue(sortedSequences[matches[i]+1], i);		//adds intex of sequence to respective group (offset by 1 is necessary since group 0 is reserved for unidentified sequences).
	}
	resizeGroups(sortedSequences);							//resizes the groups to take minimal space.
	return sortedSequences;
}

/*!
\brief Function for performing all demultiplexing operations.
\param seqs StringSet of sequences.
\param ids StringSet of IDs.
\param barcodes StringSet of barcodes.
\param approx search Boolean indication whether approximate search (true) or exact search (false) shall be applied.
\param hardClip Boolean indication wheter the first Bases of each sequence shall be clipped without considering the barcodes (true).
\return A two-dimensional vector of integers representing the indices of the sequences.
\remark The 0-th group (i.e. 0-th collum of the vector) contains the unidentified sequences.
\remark All following groups hold the sequences associated with the i-1-th barcode.
\remark On Problems with the given barcodes the vector returned holds a -1 at position [0][0];
*/
template<typename TSeqs, typename TIds, typename TBarcodes>
std::vector<std::vector<int> > DoAll(TSeqs& seqs, TIds& ids, TBarcodes& barcodes, bool approxSearch, bool hardClip)
{
	if(!check(seqs, ids, barcodes))
	{
		std::vector<std::vector<int> > sortedSequences;
		resize(sortedSequences, 1);
		appendValue(sortedSequences[0], -1);
		return sortedSequences;
	}
	TSeqs prefices = getPrefix(seqs, length(barcodes[0]));
	std::vector<int> matches;
	if(!approxSearch)		//if exact search is selected
	{
		Index<TBarcodes, IndexEsa<> > indexSet(barcodes);
		Finder<Index<TBarcodes, IndexEsa<> > > esaFinder(indexSet);
		matches = findAllExactIndex(prefices, esaFinder);
	}
	else					//if approximate search is selected
	{
		std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);
		matches = findAllApprox(seqs, patterns);
	}
	if(!hardClip)		//clip barcodes according to selected method
	{
		clipBarcodes(seqs, matches, length(barcodes[0]));
	}
	else
	{
		clipBarcodes(seqs, length(barcodes[0]));
	}
	std::vector<std::vector<int> > sortedSequences = group(seqs, matches, barcodes);
	return sortedSequences;
}

/*!
\brief Overload of ::DoAll(TSeqs& seqs, TIds& ids, TBarcodes& barcodes, bool approxSearch, bool hardClip) for paired-end data.
\param seqs StringSet of forward reads.
\param seqsRev StringSet of backward reads.
\param ids StringSet of IDs.
\param barcodes StringSet of barcodes.
\param approx search Boolean indication whether approximate search (true) or exact search (false) shall be applied.
\param hardClip Boolean indication wheter the first Bases of each sequence shall be clipped without considering the barcodes (true).
\return A two-dimensional vector of integers representing the indices of the sequences.
\remark The 0-th group (i.e. 0-th collum of the vector) contains the unidentified sequences.
\remark All following groups hold the sequences associated with the i-1-th barcode.
\remark On Problems with the given barcodes the vector returned holds a -1 at position [0][0];
*/
template<typename TSeqs, typename TIds, typename TBarcodes>
std::vector<std::vector<int> > DoAll(TSeqs& seqs, TSeqs& seqsRev, TIds& ids, TBarcodes& barcodes, bool approxSearch, bool hardClip)
{
	if(!check(seqs, seqsRev, ids, barcodes))
	{
		std::vector<std::vector<int> > sortedSequences;
		resize(sortedSequences, 1);
		appendValue(sortedSequences[0], -1);
		return sortedSequences;
	}
	TSeqs prefices = getPrefix(seqs, length(barcodes[0]));
	std::vector<int> matches;
	if(!approxSearch)		//if exact search is selected
	{
		Index<TBarcodes, IndexEsa<> > indexSet(barcodes);
		Finder<Index<TBarcodes, IndexEsa<> > >esaFinder(indexSet);
		matches = findAllExactIndex(prefices, esaFinder);
	}
	else					//if approximate search is selected
	{
		std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);
		matches = findAllApprox(seqs, patterns);
	}
	if(!hardClip)		//clip barcodes according to selected method
	{
		clipBarcodes(seqs, matches, length(barcodes[0]));
	}
	else
	{
		clipBarcodes(seqs, length(barcodes[0]));
	}
	std::vector<std::vector<int> > sortedSequences = group(seqs, matches, barcodes);
	return sortedSequences;
}

/*!
\brief Overload of ::DoAll(TSeqs& seqs, TIds& ids, TQuals& quals, TBarcodes& barcodes, bool approxSearch, bool hardClip) for multiplex barcode single-end or paired-end data.
\param seqs StringSet of forward reads.
\param ids StringSet of IDs.
\param multiplex StringSet of multiplex barcodes.
\param barcodes StringSet of barcodes.
\param approx search Boolean indication whether approximate search (true) or exact search (false) shall be applied.
\return A two-dimensional vector of integers representing the indices of the sequences.
\remark The 0-th group (i.e. 0-th collum of the vector) contains the unidentified sequences.
\remark All following groups hold the sequences associated with the i-1-th barcode.
*/
template<typename TSeqs, typename TIds, typename TBarcodes>
std::vector<std::vector<int> > DoAll(TSeqs& seqs, TIds& ids, TBarcodes& multiplex, TBarcodes& barcodes, bool approxSearch)
{
	std::vector<int> matches;
	if(!approxSearch)		//if exact search is selected
	{
		Index<TBarcodes, IndexEsa<> > indexSet(barcodes);
		Finder<Index<TBarcodes, IndexEsa<> > >esaFinder(indexSet);
		matches = findAllExactIndex(multiplex, esaFinder);
	}
	else					//if approximate search is selected
	{
		std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);
		matches = findAllApprox(multiplex, patterns);
	}
	std::vector<std::vector<int> > sortedSequences = group(seqs, matches, barcodes);
	return sortedSequences;
}
/*!
\brief Function for creating Vectors of StringSets from the vector produced by ::DoAll, oder ::group.
\param seqs StringSet of forward reads.
\param seqsRev StringSet backward reads.
\param ids StringSet of IDs.
\param groups two-dimensional vector of integers representing the indices of the sequences. Produced by ::DoAll, oder ::group.
\param gSeq Vector of StringSets to be filled with the forward reads.
\param gSeqRev Vector of StringSets to be filled with the backward reads.
\param gIds Vector of StringSets to be filled with the IDs.
\remark the given StringSet are deletet.
*/

template<typename TSeqs, typename TIds, typename TGroups>
void buildSets(TSeqs& seqs, TSeqs& seqsRev, TIds& ids, const std::vector<std::vector<int> >& groups, std::vector<TSeqs>& gSeqs, std::vector<TSeqs>& gSeqsRev, std::vector<TIds>& gIds)
{
	unsigned k = 0;
	for(unsigned i = 0; i < length(groups); ++i)
	{
		for(unsigned j = 0; j < length(groups[i]); ++j)
		{
			appendValue(gSeqs[k], seqs[group[i][j]]);
			appendValue(gSeqsRev[k], seqsRev[group[i][j]]);
			appendValue(gIds[k], ids[group[i][j]]);
		}
		if(length(groups[i]) != 0)
		{
			++k;
		}
	}
	clear(seqs);
	resize(seqs, 0);
	clear(seqsRev);
	resize(seqsRev, 0);
	clear(ids);
	resize(ids, 0);
}

/*!
\brief Overload of ::buildSets(TSeqs seqs, TSeqsRev, TIds ids, TGroups groups, std::vector<TSeqs> gSeqs, std::vector<TSeqs> gSeqsRev, std::vector<TIds> gIds) for single-end data.
\param seqs StringSet of forward reads.
\param ids StringSet of IDs.
\param groups two-dimensional vector of integers representing the indices of the sequences. Produced by ::DoAll, oder ::group.
\param gSeq Vector of StringSets to be filled with the forward reads.
\param gSeqRev Vector of StringSets to be filled with the backward reads.
\param gIds Vector of StringSets to be filled with the IDs.
\remark the given StringSet are deletet.
*/
template<typename TSeqs, typename TIds>
void buildSets(TSeqs& seqs, TIds& ids, const std::vector<std::vector<int> >& groups, std::vector<TSeqs>& gSeqs, std::vector<TIds>& gIds)
{
	resize(gSeqs, length(groups));
	resize(gIds, length(groups));
	unsigned k = 0;
	for(unsigned i = 0; i < length(groups); ++i)
	{
		for(unsigned j = 0; j < length(groups[i]); ++j)
		{
			appendValue(gSeqs[k], seqs[groups[i][j]]);
			appendValue(gIds[k],ids[groups[i][j]]);
		}
		if(length(groups[i]) != 0)
		{
			++k;
		}
	}
	clear(seqs);
	resize(seqs, 0);
	clear(ids);
	resize(ids, 0);
}


template<typename TgSeqs, typename TgIds, typename TBarcodes>
int writeGroups(TgSeqs& gSeqs, TgSeqs& gSeqsRev, TgIds& gIds, TBarcodes& barcodes, std::vector<std::vector<int> >& groups)
{
	if(writeGroups(gSeqs, barcodes, groups) != 0);		//Writing of forward reads
	{
		return 1;
	}
	unsigned j = 0;
	for(unsigned i = 0; i < length(groups); ++i)
	{
		if(lenght(groups[i]) != 0)
		{
			String<char> outputPath = "D:/seqan/out/";
			appendValue(outputPath, barcodes[i]);			//TODO: append to desired folder
			appendValue(outputPath, "_REV.fq");			
			
			SequenceStream seqStream(toCString(outputPath), SequenceStream::WRITE);
			if (!isGood(seqStream))
			{
				std::cerr << "ERROR: Could not create/open the file.\n";
				return 1;
			}
			if(writeAll(seqStream, gIds[j], gSeqs[j]) != 0)		
			{
				std::cerr << "Error while saving the sequences.\n";
				return 1;
			}
			++j;
		}
	}
	return 0;
}


template<typename TgSeqs, typename TgIds, typename TBarcodes, typename Tgroups>
int writeGroups(TgSeqs& gSeqs, TgIds& gIds, TBarcodes& barcodes, Tgroups& groups)
{
	unsigned j = 0;
	for(unsigned i = 0; i < length(groups); ++i)
	{
		if(length(groups[i]) != 0)
		{
			String<char> outputPath = "D:/seqan/out/";
			append(outputPath, barcodes[i]);			//TODO: append to desired folder
			append(outputPath, ".fq");
			SequenceStream seqStream(toCString(outputPath), SequenceStream::WRITE);
			if (!isGood(seqStream))
			{
				std::cerr << "ERROR: Could not create/open the file.\n";
				return 1;
			}
			if(writeAll(seqStream, gIds[j], gSeqs[j]) != 0)		
			{
				std::cerr << "Error while saving the sequences.\n";
				return 1;
			}
			++j;
		}
	}
	return 0;
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_
