// ==========================================================================
//                                  msplazer
// ==========================================================================
// Copyright (c) 2011, Kathrin Trappe, FU Berlin
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
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_KTRAPPE_APPS_MSPLAZER_STELLAR_ROUTINES_H_
#define SANDBOX_KTRAPPE_APPS_MSPLAZER_STELLAR_ROUTINES_H_

#define BREAKPOINT_DEBUG

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "../../../../core/apps/stellar/stellar.h"
using namespace seqan;

//Compare StellarMatches using in order of priority (1) read begin, (2) read end position, (3) chromosome begin, (4) chromosome end pos
template< typename TSequence, typename TId>
struct CompareStellarMatches{
	bool operator()
		(StellarMatch<TSequence, TId> const &  match1, StellarMatch<TSequence, TId> const & match2) const
		{
			if(match1.begin2 != match2.begin2)
				return match1.begin2 < match2.begin2;
			if(match1.end2 != match2.end2)
				return match1.end2 < match2.end2;
			if(match1.begin1 != match2.begin1)
				return match1.begin1 < match2.begin1;
			//if(match1.end1 != match2.end1)
				
			return match1.end1 < match2.end1;
		}
};

///////////////////////////////////////////////////////////////////////////////
//Wrapper functions to get Stellar matches

template <typename TSequence, typename TMatches>
void _getStellarMatches(StringSet<TSequence> & queries, StringSet<TSequence> & databases, 
		StringSet<CharString> & databaseIDs, StellarOptions & stellarOptions, TMatches & stellarMatches){

	//Finder
	typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
	//	typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;

	//Using Stellars structure for queries
	typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
	//typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape> > TQGramIndex;
	TQGramIndex qgramIndex(queries);
	//stellarOptions.qGram = 4;
	resize(indexShape(qgramIndex), stellarOptions.qGram);
	//              cargo(qgramIndex).abundanceCut = stellarOptions.qgramAbundanceCut;


	//pattern
	/*
	//Using FragmenStore
	typedef Index<StringSet<TSequence, Owner<ConcatDirect<> > >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndexFrSt;
	TQGramIndexFrSt qgramIndexFrSt(fragments.readSeqStore);
	stellarOptions.qGram = 4;
	resize(indexShape(qgramIndexFrSt), stellarOptions.qGram);
	typedef Pattern<TQGramIndexFrSt, Swift<SwiftLocal> > TPatternFrSt;
	TPatternFrSt swiftPatternFrSt(qgramIndexFrSt);
	std::cout << "FrStPattern: " << value(host(needle(swiftPatternFrSt)),0) << std::endl;
	// Construct index
	std::cout << "Constructing index..." << std::endl;
	indexRequire(qgramIndexFrSt, QGramSADir());
	std::cout << std::endl;
	*/

	//Using Stellars structure for queries
	typedef Pattern<TQGramIndex, Swift<SwiftLocal> > TPattern;
	TPattern swiftPattern(qgramIndex);
	//std::cout << "Pattern: " << value(host(needle(swiftPattern)),0) << std::endl;

	// Construct index
	std::cout << "Constructing index..." << std::endl;
	indexRequire(qgramIndex, QGramSADir());
	std::cout << std::endl;

	//Call Stellar for each database sequence
	double start = sysTime();
	for(unsigned i = 0; i < length(databases); ++i){
		//Using long stellar() to calculate stellarMatches on + strand
		if(stellarOptions.forward){
			TFinder swiftFinder(databases[i], stellarOptions.minRepeatLength, stellarOptions.maxRepeatPeriod);
			stellar(swiftFinder, swiftPattern, stellarOptions.epsilon, stellarOptions.minLength, stellarOptions.xDrop,
					stellarOptions.disableThresh, stellarOptions.compactThresh, stellarOptions.numMatches,
					stellarOptions.verbose, databaseIDs[i], true, stellarMatches, AllLocal());
		}

		// - strand
		if(stellarOptions.reverse){
			reverseComplement(databases[i]);
			TFinder revSwiftFinder(databases[i], stellarOptions.minRepeatLength, stellarOptions.maxRepeatPeriod);
			stellar(revSwiftFinder, swiftPattern, stellarOptions.epsilon, stellarOptions.minLength, stellarOptions.xDrop,
					stellarOptions.disableThresh, stellarOptions.compactThresh, stellarOptions.numMatches,
					stellarOptions.verbose, databaseIDs[i], false, stellarMatches, AllLocal());
			reverseComplement(databases[i]);
		}
	}
	std::cout << "TIME stellar " << (sysTime() - start) << "s" << std::endl;
}

//Alignment Score using _analyzeAlignment(row0, row1, alignLen, matchNumber) from stellar_output.h
//matchNumber and alignLen reset within _analyzeAlignment
template < typename TSequence, typename TId, typename TValue >
void _getScore(StellarMatch<TSequence, TId> & match, TValue & alignLen, TValue & matchNumber, TValue & alignDistance){
	
	if(!match.orientation)
		reverseComplement(infix(source(match.row1),match.begin1,match.end1));

	_analyzeAlignment(match.row1, match.row2, alignLen, matchNumber);
	alignDistance = alignLen - matchNumber;
	if(!match.orientation)
		reverseComplement(infix(source(match.row1),match.begin1,match.end1));

}

//Transforms coordinates of a stellar match onto the other strand
template <typename TMatch>
void _transformCoordinates(TMatch & match){
	
	//setClippedBeginPosition(match.row1, length(source(match.row1)) - match.end1);
	//setClippedEndPosition(match.row1, length(source(match.row1)) - match.begin1);
	//std::cerr << "Chr length: " << length(source(match.row1)) << std::endl;
	match.row1.clipped_source_begin = length(source(match.row1)) - match.end1;
	match.row1.clipped_source_end = length(source(match.row1)) - match.begin1;
	match.begin1 = clippedBeginPosition(match.row1);
	match.end1 = clippedEndPosition(match.row1);
}

//Computes distance score for each match and stores is in distanceScores
//Note: Matches are being sorted within this function
//Note also: StellarMatches of the reverse Strand are being modified, in the sense that they correspond to the right positions within the forward strand
template < typename TScoreAlloc, typename TSequence, typename TId>
void _getMatchDistanceScore(
		StringSet <QueryMatches <StellarMatch <TSequence, TId> > > & stellarMatches, 
		String<TScoreAlloc> & distanceScores, bool internalStellarCall){

	typedef StellarMatch<TSequence, TId> TMatch;
	typedef typename Size<typename TMatch::TAlign>::Type TSize;
	typedef typename Iterator<String<TMatch> >::Type TIterator;


	//Calculating distance score of each match using Stellars _analayzeAlignment function
	int matchNumber, alignLen, alignDistance;

	for(TSize i = 0; i < length(stellarMatches); ++i){

		TScoreAlloc & matchDistanceScores = distanceScores[i];
		resize(matchDistanceScores, length(stellarMatches[i].matches));
		//Sorting matches according to query begin position (begin2)
		std::sort(begin(stellarMatches[i].matches), end(stellarMatches[i].matches), CompareStellarMatches<TSequence, TId>());
		TIterator itStellarMatches = begin(stellarMatches[i].matches);
		TIterator itEndStellarMatches = end(stellarMatches[i].matches);
		TSize matchIndex = 0;	
		for(;itStellarMatches < itEndStellarMatches;goNext(itStellarMatches)){
			if(internalStellarCall && !(*itStellarMatches).orientation){
				//coordinate transformation (setClippedPosition does not work here)
				_transformCoordinates(*itStellarMatches);
			}
			//Compute edit distance score
			_getScore(*itStellarMatches, alignLen, matchNumber, alignDistance);
			matchDistanceScores[matchIndex] = alignDistance;
			
			++matchIndex;
		}
	}
}

template <typename TMatch, typename TPos>
void _trimMatchBegin(TMatch & stMatch, TPos const & splitPos, TPos const & projSplitPos){

	setClippedBeginPosition(stMatch.row2, projSplitPos);
	if(stMatch.orientation)
		setClippedBeginPosition(stMatch.row1, splitPos);
	else{

		TPos diff = splitPos - clippedBeginPosition(stMatch.row1);
		_transformCoordinates(stMatch);
		setClippedBeginPosition(stMatch.row1, stMatch.begin1 + diff);
		stMatch.begin1 = clippedBeginPosition(stMatch.row1);
		_transformCoordinates(stMatch);
	}
	stMatch.begin1 = clippedBeginPosition(stMatch.row1);
	stMatch.begin2 = clippedBeginPosition(stMatch.row2);
}

template <typename TMatch, typename TPos>
void _trimMatchEnd(TMatch & stMatch, TPos & splitPos, TPos & projSplitPos){

	setClippedEndPosition(stMatch.row2, projSplitPos);
	if(stMatch.orientation)
		setClippedEndPosition(stMatch.row1, splitPos);
	else{

		TPos diff = clippedEndPosition(stMatch.row1) - splitPos;
		_transformCoordinates(stMatch);
		setClippedEndPosition(stMatch.row1, stMatch.end1 - diff);
		stMatch.end1 = clippedEndPosition(stMatch.row1);
		_transformCoordinates(stMatch);

	}

	stMatch.end1 = clippedEndPosition(stMatch.row1);
	stMatch.end2 = clippedEndPosition(stMatch.row2);
}


//TODO Restructure
template <typename TSequence, typename TId, typename TBreakpoint>
void _getStellarIndel(StellarMatch<TSequence, TId> & match,
		String<TBreakpoint> & globalStellarIndels,
		TId const & queryId,
		TSequence & query){

		//std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	typedef typename Infix<TSequence>::Type TInfix;
	typedef typename TBreakpoint::TPos TPos;
	//typedef typename StellarMatch<TSequence, TId>::TRow TRow;
	typedef Align<TInfix> TAlign;
	typedef typename Row<TAlign>::Type TRow;

	bool gapOpen = true;
	TPos startSeqPos, endSeqPos, indelStart, pos;
	indelStart = 0;

	//Match refinement
	int matchNumber, alignLen, score;
	Score<int> scoreType(1,-2,-2,-5);

    /*	
	std::cerr << "Trying to build infices: " << length(match.row1) << " begin1: " << match.begin1 << " clippedBegin: " << clippedBeginPosition(match.row1) <<
		" end1: " << match.end1 << " clippedEnd: " << clippedEndPosition(match.row1) << std::endl;
    */	
	TInfix seq1, seq2;
	seq1 = infix(source(match.row1), match.begin1, match.end1);

	seq2 = infix(query, match.begin2, match.end2);

	if(!match.orientation)
		reverseComplement(seq1);
	
	TAlign align;
	resize(rows(align), 2);
	setSource(row(align, 0), seq1);
	setSource(row(align, 1), seq2);

	StringSet<TInfix> seqs;

	appendValue(seqs, seq1);
	appendValue(seqs, seq2);

	_analyzeAlignment(match.row1, match.row2, alignLen, matchNumber);
	score = alignLen - matchNumber;

	//Refine match using banded global alignment with affine gap score
	if(score > 0){
		globalAlignment(align, seqs, scoreType, -score, score, BandedGotoh());
		//std::cout << align << std::endl;
	}

	TRow & rowDB = row(align, 0);
	TRow & rowRead = row(align, 1);

	//Insertions
	for(pos = 0; pos < std::max(endPosition(rowDB), endPosition(rowRead)); ++pos){
		//Look for gap in rowDB
		if(isGap(rowDB,pos)){
			//Keep track if it is the first gap, i.e. the beginning of the insertion
			if(gapOpen){
				gapOpen = false;
				indelStart = pos;
				//std::cerr << "indelStart: " << indelStart << " gapopen " << gapOpen << std::endl;
			}
		}else{
			//If an insertion just closed, i.e. !gapOpen, save the insertion in a breakpoint
			if(!gapOpen){
				//Get breakpoint position (note, in case of an insertion start=end)
				startSeqPos = toSourcePosition(rowDB, indelStart) + match.begin1;
				TPos readStartPos = toSourcePosition(rowRead, indelStart) + match.begin2;
				TPos readEndPos = toSourcePosition(rowRead, pos) + match.begin2;
				TBreakpoint bp(match.id, match.id, match.orientation, match.orientation, startSeqPos, startSeqPos, readStartPos, readEndPos, queryId);
				//get insertion sequence
				TPos bPos = toSourcePosition(rowRead, indelStart);
				TPos ePos = toSourcePosition(rowRead, pos);
				if(bPos > ePos)
					std::swap(bPos, ePos);
				TInfix inSeq = infix(query, bPos, ePos);
				setSVType(bp, static_cast<TId>("insertion"));
				setInsertionSeq(bp, inSeq);
				//_insertBreakpoint(globalStellarIndels, bp);
				//std::cerr << bp << std::endl;
				//reset gapOpen
				gapOpen = true;
			}
		}
	}
	//Get insertion at the end of the match
	if(!gapOpen){
		//Get breakpoint position (note, in case of an insertion start=end)
		startSeqPos = toSourcePosition(rowDB, indelStart) + match.begin1;
		TPos readStartPos = toSourcePosition(rowRead, indelStart) + match.begin2;
		TPos readEndPos = toSourcePosition(rowRead, pos) + match.begin2;
		TBreakpoint bp(match.id, match.id, match.orientation, match.orientation, startSeqPos, startSeqPos, readStartPos, readEndPos, queryId);
		//get insertion sequence
		TPos bPos = toSourcePosition(rowRead, indelStart) + match.begin2;
		TPos ePos = toSourcePosition(rowRead, pos) + match.begin2;
		if(bPos > ePos)
			std::swap(bPos, ePos);
		TInfix inSeq = infix(query, bPos, ePos);
		setSVType(bp, static_cast<TId>("insertion"));
		setInsertionSeq(bp, inSeq);
		_insertBreakpoint(globalStellarIndels, bp);
		//reset gapOpen
		gapOpen = true;
	}

	//Deletions
	for(pos = 0; pos < std::max(endPosition(rowDB), endPosition(rowRead)); ++pos){
		
		//Look for gap in rowRead
		if(isGap(rowRead,pos)){
			if(gapOpen){
				gapOpen = false;
				indelStart = pos;
			}
		}else{
			if(!gapOpen){
				//Get breakpoint positions
				startSeqPos = toSourcePosition(rowDB, indelStart) + match.begin1;
				endSeqPos = toSourcePosition(rowDB, pos) + match.begin1;
				if(startSeqPos > endSeqPos)
					std::swap(startSeqPos, endSeqPos);
				TPos readStartPos = toSourcePosition(rowRead, indelStart) + match.begin2;
				TPos readEndPos = toSourcePosition(rowRead, pos) + match.begin2;
				TBreakpoint bp(match.id, match.id, match.orientation, match.orientation, startSeqPos, endSeqPos, readStartPos, readEndPos, queryId);
				setSVType(bp, static_cast<TId>("deletion"));
				_insertBreakpoint(globalStellarIndels, bp);
				gapOpen = true;
			}	
		}
	}
	//Get deletion at the end of the match
	if(!gapOpen){
		startSeqPos = toSourcePosition(rowDB, indelStart) + match.begin1;
		endSeqPos = toSourcePosition(rowDB, pos) + match.begin1;
		if(startSeqPos > endSeqPos)
			std::swap(startSeqPos, endSeqPos);
		TPos readStartPos = toSourcePosition(rowRead, indelStart) + match.begin2;
		TPos readEndPos = toSourcePosition(rowRead, pos) + match.begin2;
		TBreakpoint bp(match.id, match.id, match.orientation, match.orientation, startSeqPos, endSeqPos, readStartPos, readEndPos, queryId);
		setSVType(bp, static_cast<TId>("deletion"));
		_insertBreakpoint(globalStellarIndels, bp);
		gapOpen = true;
	}
	
	if(!match.orientation)
		reverseComplement(seq1);
}

template<typename TId>
bool _checkUniqueId(TId const & sId, TId const & id, StringSet<TId> & ids, StringSet<TId> & sQueryIds){

	bool unique = true;
	for(unsigned j = 0; j < length(sQueryIds); ++j){
		if(sId == sQueryIds[j]){
			std::cout << "Found nonunique sequence ID!" << std::endl;
			std::cout << ids[j] << std::endl;
			std::cout << id << std::endl;
			std::cout << "###########################" << std::endl;
			unique = false;
		}
	}
	return unique;
}

///////////////////////////////////////////////////////////////////////////////
//Functions taken from Stellar code for writing Stellar parameters and read files

///////////////////////////////////////////////////////////////////////////////
// Imports sequences from a file, 
//  stores them in the StringSet seqs and their identifiers in the StringSet ids
template<typename TSequence, typename TId>
inline bool
_importSequences(CharString const & fileName,
		CharString const & name,
		StringSet<TSequence> & seqs,
		StringSet<TId> & ids) {
	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY)) {
		std::cerr << "Failed to open " << name << " file." << std::endl;
		return false;
	}
	StringSet<TId> sQueryIds;

	AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

	unsigned seqCount = length(multiSeqFile);
	reserve(seqs, seqCount, Exact());
	reserve(ids, seqCount, Exact());

	TSequence seq;
	TId id;
	TId sId;
	unsigned counter = 0;
	for(unsigned i = 0; i < seqCount; ++i) {
		assignSeq(seq, multiSeqFile[i], format);
		assignSeqId(id, multiSeqFile[i], format);
		appendValue(seqs, seq, Generous());
		appendValue(ids, id, Generous());
		
		_getShortId(id, sId);
		if(!_checkUniqueId(sId, id, ids, sQueryIds))
			++counter;
		appendValue(sQueryIds, sId);
	}

	std::cout << "Loaded " << seqCount << " " << name << " sequence" << ((seqCount>1)?"s.":".") << std::endl;
	if(counter > 0)
    	std::cout << "Found " << counter << " nonunique sequence IDs" << std::endl;
	return true;
}


///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and writes them to std::cout
void _writeCalculatedParams(StellarOptions & options) {
//IOREV _notio_
	int errMinLen = (int) floor(options.epsilon * options.minLength);
	int n = (int) ceil((errMinLen + 1) / options.epsilon);
	int errN = (int) floor(options.epsilon * n);
	unsigned smin = (unsigned) _min(ceil((double)(options.minLength-errMinLen)/(errMinLen+1)),
		                            ceil((double)(n-errN)/(errN+1)));

	std::cout << "Calculated parameters:" << std::endl;
	if (options.qGram == (unsigned)-1) {
		options.qGram = (unsigned)_min(smin, 32u);
		std::cout << "  k-mer length: " << options.qGram << std::endl;
	}

	int threshold = (int) _max(1, (int) _min((n + 1) - options.qGram * (errN + 1),
											 (options.minLength + 1) - options.qGram * (errMinLen + 1)));
	int overlap = (int) floor((2 * threshold + options.qGram - 3) / (1 / options.epsilon - options.qGram));
	int distanceCut = (threshold - 1) + options.qGram * overlap + options.qGram;
	int logDelta = _max(4, (int) ceil(log((double)overlap + 1) / log(2.0)));
	int delta = 1 << logDelta;

	std::cout << "  s^min       : " << smin << std::endl;
	std::cout << "  threshold   : " << threshold << std::endl;
	std::cout << "  distance cut: " << distanceCut << std::endl;
	std::cout << "  delta       : " << delta << std::endl;
	std::cout << "  overlap     : " << overlap << std::endl;
	std::cout << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
// Writes user specified parameters from options object to std::cout
template<typename TOptions>
void
_writeSpecifiedParams(TOptions & options) {
//IOREV _notio_
	// Output user specified parameters
	std::cout << "User specified parameters:" << std::endl;
	std::cout << "  minimal match length             : " << options.minLength << std::endl;
	std::cout << "  maximal error rate (epsilon)     : " << options.epsilon << std::endl;
	std::cout << "  maximal x-drop                   : " << options.xDrop << std::endl;
	if (options.qGram != (unsigned)-1)
		std::cout << "  k-mer (q-gram) length            : " << options.qGram << std::endl;
	std::cout << "  search forward strand            : " << ((options.forward)?"yes":"no") << std::endl;
	std::cout << "  search reverse complement        : " << ((options.reverse)?"yes":"no") << std::endl;
	std::cout << std::endl;

	std::cout << "  verification strategy            : " << options.fastOption << std::endl;
	if (options.disableThresh != (unsigned)-1) {
		std::cout << "  disable queries with more than   : " << options.disableThresh << " matches" << std::endl;
	}
	std::cout << "  maximal number of matches        : " << options.numMatches << std::endl;
	std::cout << "  duplicate removal every          : " << options.compactThresh << std::endl;
	if (options.maxRepeatPeriod != 1 || options.minRepeatLength != 1000) {
		std::cout << "  max low complexity repeat period : " << options.maxRepeatPeriod << std::endl;
		std::cout << "  min low complexity repeat length : " << options.minRepeatLength << std::endl;
	}
	if (options.qgramAbundanceCut != 1) {
		std::cout << "  q-gram abundance cut ratio       : " << options.qgramAbundanceCut << std::endl;
	}
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes file name from options object to std::cout
template<typename TOptions>
void
_writeFileNames(TOptions & options) {
//IOREV _notio_
	std::cout << "Database file   : " << options.databaseFile << std::endl;
	std::cout << "Query file      : " << options.queryFile << std::endl;
	std::cout << "Output file     : " << options.outputFile << std::endl;
	std::cout << "Output format   : " << options.outputFormat << std::endl;
	if (options.disableThresh != (unsigned)-1) {
		std::cout << "Disabled queries: " << options.disabledQueriesFile << std::endl;
	}
	std::cout << std::endl;
}
#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_STELLAR_ROUTINES_H_

