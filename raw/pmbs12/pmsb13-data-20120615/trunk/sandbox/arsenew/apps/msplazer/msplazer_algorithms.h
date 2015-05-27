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

#ifndef SANDBOX_KTRAPPE_APPS_MSPLAZER_MSPLAZER_ALGORITHMS_H_
#define SANDBOX_KTRAPPE_APPS_MSPLAZER_MSPLAZER_ALGORITHMS_H_

#include "msplazer.h"
using namespace seqan;


//Check for match overlap
template <typename TPos>
inline bool _checkMatchOverlap(TPos const & m1Begin, TPos const & m1End, TPos const & m2Begin, TPos const & m2End){

	//if overlap: begin position of snd match is smaller than end position of fst match
	if(m2Begin < m1End && m1End < m2End && m1Begin < m2Begin)
		return true;
	return false;

}

//Check for match similarity: The overlapping part has to be smaller than a specified percentage of the each match length
template < typename TPos >
inline bool _checkMatchSim(TPos const & m1Begin, TPos const & m1End, TPos const & m2Begin,
		TPos const & m2End, MSplazerOptions const & msplazerOptions){

	double overlapPartLength = m1End - m2Begin;
	//Catch special case of required overlap=0
	if(msplazerOptions.simThresh == (double) 0.0)
		return (overlapPartLength == (double) 0.0);
	double match1Length = m1End - m1Begin;
	double match2Length = m2End - m2Begin;
	//check if overlapping percent of each match is lower than the allowed percent threshold
	if(((double) overlapPartLength/match1Length < msplazerOptions.simThresh) && ((double) overlapPartLength/match2Length < msplazerOptions.simThresh))
		return true;
	return false;

}
//Check distance between matches: The maximal allowed distance between two matches has to be smaller than a specified threshold value
//template < typename TSequence, typename TId >
template <typename TPos >
inline bool _checkMatchDist(TPos const & begin, TPos const & end, MSplazerOptions const & msplazerOptions){
	
	//Check if distance between matches is smaller than the distance threshold, assumes that begin > end
	if((int) (begin - end) < (msplazerOptions.gapThresh + 1))
		return true;
	return false;
}

//Check matches for same database
template < typename TSequence, typename TId >
inline bool _checkDBIds(StellarMatch <TSequence, TId> const & match1, StellarMatch <TSequence, TId> const & match2){

	return (match1.id == match2.id);
}

//Check matches for same strand
template < typename TSequence, typename TId >
inline bool _checkMatchStrands(StellarMatch <TSequence, TId> const & match1, StellarMatch <TSequence, TId> const & match2){

	return (match1.orientation == match2.orientation);
}

//Check order in matching strand: Function assumes that both matches are on the same strand and in the same genome and that match1.begin2 < match2.begin2 (ordered matches according to read)
template < typename TSequence, typename TId >
inline bool _checkMatchOrderInDB(StellarMatch <TSequence, TId> const & match1, StellarMatch <TSequence, TId> const & match2){

	return (match1.begin1 < match2.begin1);
}

//Check match compatibility: Checks compatibility property of two stellar matches: Do they overlap (in a specified way)? If not, are they still close enough? Returns false, if not
template <typename TPos>
inline bool _checkMatchComp(TPos const & m1Begin, TPos const & m1End, TPos const & m2Begin, TPos const & m2End,
		bool & doBP, bool & insertEdge, MSplazerOptions const & msplazerOptions){

	doBP = false;
	insertEdge = false;
	//Check for true overlap
	if(_checkMatchOverlap(m1Begin, m1End, m2Begin, m2End)){
		//Check for similarity, i.e. the percentage of the overlapping part to the match length
		doBP = _checkMatchSim(m1Begin, m1End, m2Begin, m2End, msplazerOptions);
		insertEdge = doBP;
		return true;
	}
	//If not overlapping correctly, check if there is a gap and then the gap length (matchDist)
	if(m1End < m2Begin){
		insertEdge = _checkMatchDist(m2Begin, m1End, msplazerOptions);
		//If gap length is small enough, next match might still be ok. If this gap is too big, the next one will be as well
		return insertEdge;
	}
	return true;
}

//Intitialisation of graph structure for combinable StellarMatches of a read
template < typename TSequence, typename TId, typename TGraph, typename TScoreAlloc, typename TVertexDescriptor, typename TBreakpointMap >
void _initialiseGraph(QueryMatches< StellarMatch <TSequence, TId > > & queryMatches, TGraph & graph, TScoreAlloc & matchDistanceScores, 
		TVertexDescriptor & startVertex, TVertexDescriptor & endVertex, TBreakpointMap & queryBreakpoints, MSplazerOptions const & msplazerOptions){

	//std::cerr << " Initialising graph structure " << std::endl;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<String<StellarMatch<TSequence, TId > > >::Type TIterator;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	TIterator itStellarMatches = begin(queryMatches.matches);
	TIterator itEndStellarMatches = end(queryMatches.matches);
	//The default vertex descriptor is an integer so if inserted in the same order the vertex descriptor value is the same as the position of the corresponding
	//vertex within the QueryMatches --> since we can easily iterate through the QueryMatches and use the iterator we wont keep track of the vertex descriptors
	for(;itStellarMatches < itEndStellarMatches;goNext(itStellarMatches)) {
		addVertex(graph);
	}	

	//std::cerr << " Created graph " << std::endl;
	//Add start and end to graph and property map
	startVertex = addVertex(graph);
	endVertex = addVertex(graph);

	int cargo = 0;
	resize(queryBreakpoints.slotLookupTable, 2*length(queryMatches.matches));
	//Adding edges to start vertex
	for(unsigned i = 0; i < length(queryMatches.matches); ++i){

		cargo = static_cast<int>(queryMatches.matches[i].begin2) + matchDistanceScores[i];
		if(cargo < (msplazerOptions.initGapThresh + 1)){
			TEdgeDescriptor edge = addEdge(graph, startVertex, i, cargo);//10u);//value(itV));
			resizeEdgeMap(graph, queryBreakpoints.slotLookupTable);
			assignProperty(queryBreakpoints, edge);
			//std::cerr << " Edge from start to " << i << " with cargo: " << cargo << std::endl;
		}
		cargo = static_cast<int>(length(source(queryMatches.matches[i].row2))) - static_cast<int>(queryMatches.matches[i].end2);
		if(cargo < (msplazerOptions.initGapThresh + 1)){
			TEdgeDescriptor edge = addEdge(graph, i, endVertex, cargo);
			resizeEdgeMap(graph, queryBreakpoints.slotLookupTable);
			assignProperty(queryBreakpoints, edge);
			//std::cerr << " Edge from: " << i << " to end with cargo: " << cargo << std::endl;
		}
	}
	//::std::cout << graph << ::std::endl;
}

//Match Chaining for one query: Inserts edges between compatible matches and determines their breakpoint
template < typename TSequence, typename TId, typename TGraph, typename TScoreAlloc, typename TBreakpointMap >
void _chainMatches(QueryMatches<StellarMatch<TSequence, TId > > & queryMatches,
		TId const & queryId,
		TSequence & query,
		TGraph & graph,
		TScoreAlloc & matchDistanceScores, 
		TBreakpointMap & queryBreakpoints,
		MSplazerOptions const & msplazerOptions,
		String<unsigned> & insertionCount,
		unsigned long & edgeCount,
		unsigned & noneCount){

	//For testing
	bool overlapBP = false;

	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Align<TSequence, ArrayGaps> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<TSequence>::Type	TPos;
	typedef Breakpoint<TSequence,TId> TBreakpoint;
	typedef typename Infix<TSequence>::Type TInfix;
	
	//Output values for compatibility check: do breakpoint evaluation, insert edge into graph, found gap --> stop iterating if gap is too big (overlap in read)
	bool doBP, insertEdge;
	//Penalties
	int diffDBPen, diffStrandPen, diffOrderPen;
	//Terminating condition for taking the next snd match for comparison: takeNextMatch == false meaning this and all the next matches are too far away
	bool takeNextMatch = true;
	//Breakpoint parameters
	bool trimMatches = false;
	// edit distance
	Score<int> scoreType(0,-1,-1,-1);
	int splitPos = 0;
	
	//Cargo on edges
	int cargo = 0;
	//Sequence number in match1 and match2 of query/read sequence (should be equal for all (Stellar!)matches)
	int readRowNum = 1;
	int score = 0;
	//loop over all query matches
	//std::cerr << "In chainQueryMatches length(queryMatches.matches): " << length(queryMatches.matches) << std::endl;
	for(unsigned m1 = 0; m1 < (length(queryMatches.matches) - 1); ++m1){
		
		StellarMatch< TSequence, TId > & stMatch1 = queryMatches.matches[m1];
		TPos m1Begin = stMatch1.begin2;
		TPos m1End = stMatch1.end2;
		takeNextMatch = true;
		//loop over all query matches that are supposed to be compatible
		for(unsigned m2 = m1+1; takeNextMatch && (m2 < length(queryMatches.matches)); ++m2){

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Compatibility check
			StellarMatch< TSequence, TId > & stMatch2 = queryMatches.matches[m2];
			TPos m2Begin = stMatch2.begin2;
			TPos m2End = stMatch2.end2;
			//Returns false, if next match is def. not compatible anymore
			takeNextMatch = _checkMatchComp(m1Begin, m1End, m2Begin, m2End, doBP, insertEdge, msplazerOptions);
		
			//match is compatible
			if(insertEdge){

				overlapBP = false;
				
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//Penalties
				//Different reference sequence penalty
				diffDBPen = _checkDBIds(stMatch1, stMatch2) ? 0 : msplazerOptions.diffDBPen;
				//Different orientation penalty
				diffStrandPen = (diffDBPen > 0 || _checkMatchStrands(stMatch1, stMatch2)) ? 0 : msplazerOptions.diffStrandPen;
				//Different order in reference than in read penalty
				diffOrderPen = (diffDBPen > 0 || diffStrandPen > 0 || _checkMatchOrderInDB(stMatch1, stMatch2)) ? 0 : msplazerOptions.diffOrderPen;

				//std::cerr << "doBP: " << doBP << " insertEdge: " << insertEdge << " takeNextMatch: " << takeNextMatch << std::endl;

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//Breakpoint computation, Graph input	
				//Do breakpoint computation if flag doBP was set 'true' during match compatibility check

				//Compute breakpoint
				TPos startSeqPos, endSeqPos, readStartPos, readEndPos;
				if(doBP){
					//Create alignments from Stellar rows as input for breakpoint function
					TAlign match1, match2;
					appendValue(rows(match1), stMatch1.row1);
					appendValue(rows(match1), stMatch1.row2);
					appendValue(rows(match2), stMatch2.row1);
					appendValue(rows(match2), stMatch2.row2);

					//std::cerr << match1 << match2 << std::endl;

					//Compute breakpoint score
    	    	    score = combineMatchPair(match1, match2, readRowNum, readRowNum, !stMatch1.orientation, !stMatch2.orientation, scoreType, splitPos, trimMatches);
					//Compute cargo, reduce distance by score
					cargo = static_cast<int>(matchDistanceScores[m2]) - score + diffDBPen + diffStrandPen + diffOrderPen;
					startSeqPos = toSourcePosition(stMatch1.row1,toViewPosition(stMatch1.row2, splitPos));
    	    		endSeqPos = toSourcePosition(stMatch2.row1,toViewPosition(stMatch2.row2, splitPos));
					readStartPos = splitPos;
					readEndPos = splitPos;

					++insertionCount[0];
					overlapBP = true;
				}
				else{
					//Compute score, add gap length to distance
					cargo = static_cast<int>(matchDistanceScores[m2]) + static_cast<int>(m2Begin - m1End)  + diffDBPen + diffStrandPen + diffOrderPen;
					startSeqPos = stMatch1.end1;
    	    		endSeqPos = stMatch2.begin1;
					readStartPos = stMatch1.end2;
					readEndPos = stMatch2.begin2;

					++insertionCount[static_cast<int>(m2Begin - m1End)];
				}
				++edgeCount;

				//TBreakpoint bp(stMatch1.id, stMatch2.id, stMatch1.orientation, stMatch2.orientation, startSeqPos, endSeqPos, queryId);
				TBreakpoint bp(stMatch1.id, stMatch2.id, stMatch1.orientation, stMatch2.orientation, startSeqPos, endSeqPos, readStartPos, readEndPos, queryId, overlapBP);
				//Returns true for insertion type, get insertion infix then
				if(setSVType(bp)){
				   	if(stMatch1.end2 < stMatch2.begin2){
						//TSequence inSeq;
						TInfix inSeq;
						//get insertion sequence from matches and read sequence --> infix mit endPos(match1) and startPos(match2)
						inSeq = infix(query, stMatch1.end2, stMatch2.begin2);
						//std::cerr << "Insertion sequence: " << inSeq << std::endl;
						setInsertionSeq(bp, inSeq);
						//std::cerr << "Insertion sequence in bp: " << bp.insertionSeq << std::endl;
					}else{
						//Double overlap check (not handled jet)
						//std::cerr << "double overlap in reference and read called from read overlap" << std::endl;
						clear(bp.svtype);
						bp.svtype = "none";
						++noneCount;
					}
				}	
				//std::cerr << " Edge from: " << m1 << " to " << m2 << " with cargo: " << cargo << std::endl;
				//std::cerr << " matchDistScore: " << static_cast<int>(matchDistanceScores[m2]) << " score: " << score << " diffDBPen: " << diffDBPen << " diffStrandPen: " << diffStrandPen << std::endl;
					
				//Insert breakpoint
				TEdgeDescriptor edge = addEdge(graph, m1, m2, cargo);
				resizeEdgeMap(graph, queryBreakpoints.slotLookupTable);
				assignProperty(queryBreakpoints, edge, bp);
				//std::cerr << bp << std::endl;			
			}
			doBP = false;
			insertEdge = false;
		}
	}

}


//Match Chaining for one query: Inserts edges between compatible matches and determines their breakpoint
template < typename TSequence, typename TId, typename TGraph, typename TScoreAlloc, typename TBreakpointMap >
void _chainMatchesReference(QueryMatches<StellarMatch<TSequence, TId > > & queryMatches, 
		TId const & queryId,
		TSequence & query,
		TGraph & graph,
		TScoreAlloc & matchDistanceScores,
		TBreakpointMap & queryBreakpoints,
		MSplazerOptions const & msplazerOptions,
		String<unsigned> & insertionCount,
		unsigned long & edgeCount,
		unsigned & replacedBPCount,
		unsigned & replacedInsertBPCount,
		unsigned & noneCount){

	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Align<TSequence, ArrayGaps> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<TSequence>::Type	TPos;
	typedef Breakpoint<TSequence,TId> TBreakpoint;
	typedef StellarMatch<TSequence,TId> TMatch;
	typedef typename Infix<TSequence>::Type TInfix;
	
	//Penalties
	//int diffOrderPen, diffStrandPen, insertionPen;
	int diffStrandPen, diffOrderPen;
	//Output values for compatibility check: do breakpoint evaluation, insert edge into graph, found gap --> stop iterating if gap is too big (overlap in read)
	bool doBP, insertEdge, swap;
	//Breakpoint parameters
	bool trimMatches = false;
	// edit distance
	Score<int> scoreType(0,-1,-1,-1);
	int splitPos = 0;
	//Cargo on edges
	int cargo = 0;

	//Sequence number in match1 and match2 of query/read sequence (should be equal for all (Stellar!)matches)
	int readRowNum = 1;
	//Breakpoint score (gain of edit distance)
	int score = 0;
	//loop over all query matches
	for(unsigned m1 = 0; m1 < (length(queryMatches.matches) - 1); ++m1){
		
		TMatch *stMatch1 = & queryMatches.matches[m1];
		TPos m1Begin = (*stMatch1).begin1;
		TPos m1End = (*stMatch1).end1;
		//loop over all query matches that are supposed to be compatible
		for(unsigned m2 = m1+1; m2 < length(queryMatches.matches); ++m2){

			swap = false;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Compatibility check
			TMatch *stMatch2 = & queryMatches.matches[m2];

			//std::cerr << "stMatch1 and stMatch2: " << *stMatch1 << *stMatch2 << std::endl;
			//Check if on same chromosome
			if(_checkDBIds(*stMatch1,*stMatch2)){
				//std::cout << "Matches on same reference, doing match compatibility check" << std::endl;
				//Exchange matches once, via temp reference, then change back at the end!
				if(m1Begin > (*stMatch2).begin1){
					std::swap(stMatch1,stMatch2);
					m1Begin = (*stMatch1).begin1;
					m1End = (*stMatch1).end1;  
					swap = true;
				}

				TPos m2Begin = (*stMatch2).begin1;
				TPos m2End = (*stMatch2).end1;
				_checkMatchComp(m1Begin, m1End, m2Begin, m2End, doBP, insertEdge, msplazerOptions);

				if(insertEdge){
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//Penalties
					//Insertion penalty
					//insertionPen = msplazerOptions.insertionPen + insertionLength;
					//Different orientation of matches
					diffStrandPen = (_checkMatchStrands(*stMatch1, *stMatch2)) ? 0 : msplazerOptions.diffStrandPen;
					//Different order in reference than in read penalty
					diffOrderPen = 0;
					if(swap)
						diffOrderPen = (diffStrandPen > 0 || _checkMatchOrderInDB(*stMatch1, *stMatch2)) ? 0 : msplazerOptions.diffOrderPen;
					//std::cerr << "doBP: " << doBP << " insertEdge: " << insertEdge << std::endl;

					///////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//Breakpoint computation, Graph input	
					//Do breakpoint computation if flag doBP was set 'true' during match compatibility check

					//Compute breakpoint
					TPos startSeqPos, endSeqPos, readStartPos, readEndPos;
					if(doBP){
						//std::cout << "Doing bp computation" << std::endl;
						//Create alignments from Stellar rows as input for breakpoint function
						TAlign match1, match2;
						appendValue(rows(match1), (*stMatch1).row1);
						appendValue(rows(match1), (*stMatch1).row2);
						appendValue(rows(match2), (*stMatch2).row1);
						appendValue(rows(match2), (*stMatch2).row2);
						//std::cerr << match1 << match2;
		
						//Compute breakpoint score
        	    		score = combineMatchPair(match1, match2, 1 - readRowNum, 1 - readRowNum, !(*stMatch1).orientation, !(*stMatch2).orientation, scoreType, splitPos, trimMatches);
						//std::cerr << match1 << match2;
						//Compute cargo, reduce distance by score
						cargo = static_cast<int>(matchDistanceScores[m2]) - score + diffStrandPen + diffOrderPen;//+ insertionPen
						//startSeqPos = toSourcePosition((*stMatch1).row1,toViewPosition((*stMatch1).row1, splitPos - 1));
						startSeqPos = toSourcePosition((*stMatch1).row1,toViewPosition((*stMatch1).row1, splitPos));
	    	    		endSeqPos = startSeqPos;
						//endSeqPos = toSourcePosition((*stMatch2).row1,toViewPosition((*stMatch2).row1, splitPos));
						if(!swap){
							readStartPos = toSourcePosition((*stMatch1).row2, toViewPosition((*stMatch1).row1, splitPos));
							readEndPos = toSourcePosition((*stMatch2).row2, toViewPosition((*stMatch2).row1, splitPos));
						}else{
							readStartPos = toSourcePosition((*stMatch2).row2, toViewPosition((*stMatch2).row1, splitPos));
							readEndPos = toSourcePosition((*stMatch1).row2, toViewPosition((*stMatch1).row1, splitPos));
						}
					}
					else{
						//Compute score, add gap length to distance
						cargo = static_cast<int>(matchDistanceScores[m2]) + static_cast<int>(m2Begin - m1End) + diffStrandPen + diffOrderPen;//+ insertionPen
						if(!swap){
							startSeqPos = (*stMatch1).end1;//m1End
	    	    			endSeqPos = (*stMatch2).begin1;
							readStartPos = (*stMatch1).end2;
							readEndPos = (*stMatch2).begin2;
						}else{
							startSeqPos = (*stMatch2).end1;//m1End
    	    				endSeqPos = (*stMatch1).begin1;
							readStartPos = (*stMatch2).end2;
							readEndPos = (*stMatch1).begin2;
						}

					}
					TBreakpoint bp((*stMatch1).id, (*stMatch2).id, (*stMatch1).orientation, (*stMatch2).orientation, startSeqPos, endSeqPos, readStartPos, readEndPos, queryId);
					//std::cerr << " split pos: " << splitPos << std::endl;
					//std::cerr << "bp: " << bp << std::endl;
					
					//if(setSVType(bp) && ((*stMatch1).end2 < (*stMatch2).begin2)){
					if(setSVType(bp)){
							//TSequence inSeq;
							TInfix inSeq;
							TPos bPos, ePos;
							//get insertion sequence from matches and read sequence --> infix mit endPos(match1) and startPos(match2)
							//std::cerr << *stMatch1 << *stMatch2 << std::endl;
							if(!swap){
								//bPos = toSourcePosition((*stMatch1).row2,toViewPosition((*stMatch1).row1, splitPos));
    				    		//ePos = toSourcePosition((*stMatch2).row2,toViewPosition((*stMatch2).row1, splitPos));
								bPos = toSourcePosition((*stMatch1).row2,toViewPosition((*stMatch1).row1, splitPos));
    				    		ePos = toSourcePosition((*stMatch2).row2,toViewPosition((*stMatch2).row1, splitPos));
							}else{
								ePos = toSourcePosition((*stMatch1).row2,toViewPosition((*stMatch1).row1, splitPos));
    				    		bPos = toSourcePosition((*stMatch2).row2,toViewPosition((*stMatch2).row1, splitPos));
							}

							//std::cerr << *stMatch1 << *stMatch2 << std::endl;
							if(bPos < ePos)
								inSeq = infix(query, bPos, ePos);
							else
								inSeq = infix(query, ePos, bPos);
							//std::cerr << "splitPos: " << splitPos << " bPos: " << bPos << " ePos: " << ePos << " length(inSeq): " << length(inSeq) << std::endl;
							//std::cerr << "Insertion sequence: " << inSeq << std::endl;
							if(length(inSeq) == 0)
								setSVType(bp, static_cast<TId>("none"));
							else
								setInsertionSeq(bp, inSeq);
					}	
				
					//std::cerr << " Reference Edge from: " << m1 << " to " << m2 << " with cargo: " << cargo << std::endl;
					//std::cerr << " matchDistScore: " << static_cast<int>(matchDistanceScores[m2]) << " score: " << score << " diffStrandPen: " << diffStrandPen << std::endl;
					
					TEdgeDescriptor edge;
					edge = findEdge(graph, m1, m2);
					//Returns 0 if edge does not exist
					if(edge != 0){
						//Replace cargo and bp
						if(cargo < getCargo(edge)){
							assignCargo(edge, cargo);
							TBreakpoint & oldBp = property(queryBreakpoints, edge);
							//std::cerr << " REPLACE BREAKPOINT: " << oldBp << std::endl;
							if(!oldBp.overlapBP)
								++replacedInsertBPCount;
							if(getSVType(oldBp) == "none")
								--noneCount;
							oldBp = bp;
							//std::cerr << " NEW BREAKPOINT: " << oldBp << std::endl;
							++replacedBPCount;
						}
					}
					else{
						//Insert breakpoint
						edge = addEdge(graph, m1, m2, cargo);
						resizeEdgeMap(graph, queryBreakpoints.slotLookupTable);
						assignProperty(queryBreakpoints, edge, bp);
						//std::cerr << bp << std::endl;			
						if(doBP)
							++insertionCount[0];
						else
							++insertionCount[static_cast<int>(m2Begin - m1End)];

						++edgeCount;
					}
				}
				if(swap){
					std::swap(stMatch1,stMatch2);
					m1Begin = (*stMatch1).begin1;
					m1End = (*stMatch1).end1;  
					//std::cout << "Switching matches back" << std::endl;
				}

			}
		}
	}

}

//Chain all matches of each query
template < typename TSequence, typename TId, typename TScoreAlloc, typename TMSplazerChain >
void _chainQueryMatches(StringSet <QueryMatches< StellarMatch <TSequence, TId > > > & stellarMatches,
		String< TScoreAlloc> & distanceScores,
		String< TMSplazerChain > & queryChains,
		StringSet<TId> const & queryIds,
		StringSet<TSequence> & queries,
		MSplazerOptions const & msplazerOptions){

	//Debug and statistic variables
	unsigned replacedBPCount = 0;
	unsigned replacedInsertionBPCount = 0;
	unsigned emptyChainCount = 0;
	String<unsigned> insertionCount;
	resize(insertionCount, msplazerOptions.gapThresh + 1);
	for(unsigned i = 0; i < length(insertionCount); ++i)
		insertionCount[i] = 0;
	String<unsigned> insertionCountRef;
	resize(insertionCountRef, msplazerOptions.gapThresh + 1);
	for(unsigned i = 0; i < length(insertionCountRef); ++i)
		insertionCountRef[i] = 0;

	unsigned long edgeCount = 0;
	unsigned noneCount = 0;
	
	for(unsigned i = 0; i < length(stellarMatches); ++i){

		TScoreAlloc matchDistanceScores = distanceScores[i];
		TMSplazerChain chain(matchDistanceScores);
		//chain.queryId = queryIDs[i];
		//chain.databaseId = stellarMatches[i].

		//if(length(stellarMatches[i].matches) < 2){
		if(length(stellarMatches[i].matches) == 0){
			//Insert single match into graph, no extra edges beside from start and to end
			chain.isEmpty = true;	
			++emptyChainCount;
		}else{

			//Graph init
			_initialiseGraph(stellarMatches[i], 
					chain.graph, 
					chain.matchDistanceScores, 
					chain.startVertex, 
					chain.endVertex,
					chain.breakpoints,
					msplazerOptions);


			//Chain compatible matches
			_chainMatches(stellarMatches[i], 
					queryIds[i],
					queries[i],
					chain.graph, 
					chain.matchDistanceScores,
				    chain.breakpoints,	
					msplazerOptions,
					insertionCount,
					edgeCount,
					noneCount);
			
			_chainMatchesReference(stellarMatches[i], 
					queryIds[i],
					queries[i],
					chain.graph, 
					chain.matchDistanceScores,
				    chain.breakpoints,	
					msplazerOptions,
					insertionCountRef,
					edgeCount,
					replacedBPCount,
					replacedInsertionBPCount,
					noneCount);
			
		}
		appendValue(queryChains,chain);
	}
	/*
	std::cerr << "Number of edges: " << edgeCount << std::endl;
	std::cerr << "Empty chain count: " << emptyChainCount << std::endl;
	for(unsigned i = 0; i < length(insertionCount); ++i)
		std::cerr << "Number of edges with a " << i << " bp long insertion: " << insertionCount[i] << std::endl;
	for(unsigned i = 0; i < length(insertionCountRef); ++i)
		std::cerr << "REF Number of edges with a " << i << " bp long insertion: " << insertionCountRef[i] << std::endl;
	std::cerr << "Replaced BP Count: " << replacedBPCount << std::endl;
	std::cerr << "Replaced INSERTION BP Count: " << replacedInsertionBPCount << std::endl;
	std::cerr << "NONE count: " << noneCount << std::endl;
	*/
}

//Analyze chains in read graphs
template < typename TMSplazerChain>
void _analyzeChains(String<TMSplazerChain> & queryChains){

	InternalMap<int> weightMap;
	typedef typename TMSplazerChain::TGraph TGraph;
	typedef typename Size<TGraph>::Type TGraphSize;
	//String<TVertexDescriptor> predMap;
	//String<TGraphSize> distMap;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	//std::cerr << " length(queryChains) " << length(queryChains) << std::endl;
	

	for(unsigned i = 0; i < length(queryChains); ++i){

		if(!queryChains[i].isEmpty){


			resizeVertexMap(queryChains[i].graph, queryChains[i].distMap);
			dagShortestPath(queryChains[i].graph, queryChains[i].startVertex, weightMap, queryChains[i].predMap, queryChains[i].distMap);

		}else{
			//std::cerr << " Chain for query " << i << ", is empty!" << std::endl;//", queryID: " << queryIDs[i] << 
		}
	}
}

//Insert Breakpoint into string of breakpoints if it is not already in the set
//Returns true if breakpoint was new and has been inserted or false if breakpoint was already in the set (and just has been counted)
//..note: Two breakpoints are considered equal if belong to the same two sequences and if both positions are the same resp. and they have the same type. 
//For insertions, also the insertion length has to be the same
template <typename TBreakpoint>
bool _insertBreakpoint(String<TBreakpoint> & countedBP, TBreakpoint const & bp){
	
	//Breakpoint bp is compared to each breakpoint in the list (tempBP)
	//std::cout << "bp: " << bp << std::endl;
	for(unsigned i = 0; i < length(countedBP); ++i){
		TBreakpoint & tempBP = countedBP[i];
		//std::cout << "tempbp: " << i << " " << tempBP << std::endl;
		//Breakpoint comparison
		if(bp == tempBP){
			//std::cout << "Found same breakpoint! " << std::endl;
			//add new supporting Ids, automatically sets new support value
			appendSupportId(tempBP, bp.supportIds);
			//++tempBP.support;
			//++bp.support;
			return false;
		}
	}
	//std::cout << "Append new BP " << std::endl;
	appendValue(countedBP, bp);
	return true;
}

template <typename TBreakpoint>
void _insertBreakpoints(String<TBreakpoint> & countedBP, String<TBreakpoint> const & newBP){

	for(unsigned i = 0; i < length(newBP); ++i)
		_insertBreakpoint(countedBP, newBP[i]);
}

template <typename TMatch, typename TPos>
void _trimMatches(TMatch & fstMatch, TMatch & sndMatch, String<TPos> & splitPos){

	TPos & readEndSplitPos = splitPos[2];
	TPos & readBeginSplitPos = splitPos[3];
	TPos & refEndSplitPos = splitPos[0];
	TPos & refBeginSplitPos = splitPos[1];

	_trimMatchBegin(sndMatch, refBeginSplitPos, readBeginSplitPos);
	_trimMatchEnd(fstMatch, refEndSplitPos, readEndSplitPos);
}

template <typename TMatch, typename TPos>
void _trimMatches(String<TMatch> & matchChain, StringSet<String<TPos> > & splitPos){

	for(unsigned i = 0; i < length(splitPos); ++i){
	
		//Only trim matches that overlap in either reference or read sequence	
		if(_checkMatchOverlap(matchChain[i+1].begin1, matchChain[i+1].end1, matchChain[i].begin1, matchChain[i].end1)
				|| _checkMatchOverlap(matchChain[i+1].begin2, matchChain[i+1].end2, matchChain[i].begin2, matchChain[i].end2))
			_trimMatches(matchChain[i+1], matchChain[i], splitPos[i]);
	}
}

template <typename TBreakpoint, typename TSequence, typename TId, typename TMatch>
void _getChainIndels(String<String<TMatch> > & bestChains, String<TBreakpoint> & globalStellarIndels, TId const & queryId, TSequence & query){


	typedef typename Iterator<String<TMatch> >::Type TIterator;
	for(unsigned i = 0; i < length(bestChains); ++i){
		TIterator itStellarMatches = begin(bestChains[i]);
		TIterator itEndStellarMatches = end(bestChains[i]);
		for(;itStellarMatches < itEndStellarMatches;goNext(itStellarMatches))
			_getStellarIndel(*itStellarMatches, globalStellarIndels, queryId, query);
	}
}

template<typename TMSplazerChain>
void _findPartialChains(TMSplazerChain & queryChain){
	
	if(outDegree(queryChain.graph, queryChain.startVertex) > 0 && inDegree(queryChain.graph, queryChain.endVertex) > 0)
		return;

	if(inDegree(queryChain.graph, queryChain.endVertex) > 0){
		queryChain.isPartial = true;
		return;
	}

	if(outDegree(queryChain.graph, queryChain.startVertex) > 0){
		queryChain.isPartial = true;
		return;
	}
	return;
}

//Finding the best chain (belonging to the shortest path) and reporting the breakpoints, if any.
template < typename TBreakpoint, typename TMSplazerChain, typename TMatch>
bool _findBestChain(TMSplazerChain  & queryChain, String<TMatch> & stellarMatches,
		String<TBreakpoint> & tmpGlobalBreakpoints, MSplazerOptions const & msplazerOptions, unsigned & bcc){
		//String<TBreakpoint> & tmpGlobalBreakpoints, String<TBreakpoint> & tmpStellarIndels, unsigned & bcc){

	typedef typename TMSplazerChain::TVertexDescriptor TVertexDescriptor;
	typedef typename TMSplazerChain::TGraph TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename TBreakpoint::TPos TPos;
	typedef typename TBreakpoint::TId TId;
	bool foundBP = false;
	TVertexDescriptor spVertex1, spVertex2;
	TEdgeDescriptor edge;

	String<TMatch> tmpBestChain;
	StringSet<String<TPos> > tmpSplitPos;

	if(queryChain.isEmpty)
		return 0;
	spVertex1 = queryChain.endVertex;
	//std::cerr << "Length of stellarMatches: " << length(stellarMatches) << std::endl;
	while(spVertex1 != queryChain.startVertex){
		//Getting vertex descriptor of anchestor
		spVertex2 = getProperty(queryChain.predMap, spVertex1);
		//if vertex descriptor has max value then there is no (shortest) path from start to end.
		if(spVertex2 == maxValue<TVertexDescriptor>()){
			//queryChain.isPartial = true;
			++bcc;
			return 0;			
		}
		//Append tmpStellarMatch
		//std::cerr << "ID: " << _getId(spVertex2) << std::endl;
		if(_getId(spVertex2) != queryChain.startVertex){
			TMatch stMatch = stellarMatches[_getId(spVertex2)];
			appendValue(tmpBestChain, stMatch);
		}
		//Getting edge between both vertices
		edge = findEdge(queryChain.graph, spVertex2, spVertex1);
		// If you ever want just a copy of the breakpoints in the global set, redefine getProperty to give back getValue() instead of value(),
		// keep in mind that bp has to be defined somewhere else or it gets lost
		TBreakpoint bp;

		//Get breakpoint on edge
		if(edge != 0)
			foundBP = getProperty(queryChain.breakpoints, edge, bp);
		//else
			//std::cerr << " Edge does not exist, which makes no sense, since you are on the shortest path!" << std::endl;

		//Append breakpoint to tmpGlobalBreakpoints
		if(foundBP){
			//insert BP into tmpGlobalBreakpoints
			appendValue(tmpGlobalBreakpoints, bp);
			foundBP = false;				
			//Append splitPos
			String<TPos> splitPos;
			resize(splitPos, 4);
		
		//std::cerr << bp << std::endl;

			//Store breakpoint positions, mind. matches of different oder for translocations and reverse strand deletions
			if((getSVType(bp) == static_cast<TId>("translocation") && bp.startSeqId == bp.endSeqId) || bp.revStrandDel){
				splitPos[0] = bp.endSeqPos;
				splitPos[1] = bp.startSeqPos;
				splitPos[2] = bp.readStartPos;
			   	splitPos[3] = bp.readEndPos;
			}else{	
				splitPos[0] = bp.startSeqPos;
				splitPos[1] = bp.endSeqPos;
				splitPos[2] = bp.readStartPos;
			   	splitPos[3] = bp.readEndPos;
			}
			appendValue(tmpSplitPos, splitPos);
		}
		spVertex1 = spVertex2;
	}
	
		
	//only if chain was complete
	//trim matches with split pos
	if(msplazerOptions.simThresh > 0)
		_trimMatches(tmpBestChain, tmpSplitPos);
	
	//insert bestChain into queryChain object
	insertBestChain(queryChain, tmpBestChain);
	return 1;

}

//Finding the best chain (belonging to the shortest path) and reporting the breakpoints, if any.
template < typename TMSplazerChain, typename TBreakpoint, typename TQueryMatches, typename TSequence, typename TId>
void _findAllBestChains(String<TMSplazerChain> & queryChains,
		StringSet <TQueryMatches> & queryMatches,
		StringSet<TSequence> & queries,
		StringSet<TId> const & queryIds,
		String<TBreakpoint> & globalBreakpoints,
		String<TBreakpoint> & globalStellarIndels,
		MSplazerOptions const & msplazerOptions){

	/*
	String<unsigned> chainSizeCount;
	resize(chainSizeCount, 10);
	for(unsigned i = 0; i < length(chainSizeCount); ++i)
		chainSizeCount[i] = 0;
	*/
	typedef typename TMSplazerChain::TMatchAlloc TMatchAlloc;
	unsigned brokenChainCount = 0;
	for(unsigned i = 0; i < length(queryChains); ++i){
	
		String<TBreakpoint> tmpGlobalBreakpoints;
		//String<TBreakpoint> tmpStellarIndels;
		if(_findBestChain(queryChains[i], queryMatches[i].matches, tmpGlobalBreakpoints, msplazerOptions, brokenChainCount)){
			_insertBreakpoints(globalBreakpoints, tmpGlobalBreakpoints);
			//get small indels from matches
			_getChainIndels(queryChains[i].bestChains, globalStellarIndels, queryIds[i], queries[i]);
			/*
			for(unsigned j = 0; j < length(queryChains[i].bestChains); ++j)
				++chainSizeCount[length(queryChains[i].bestChains[j])];
		    */
		}
	}
	/*
	for(unsigned i = 0; i < length(chainSizeCount); ++ i){
		std::cout << "Number of chains with matches " << i << " : " << chainSizeCount[i] << std::endl;
	}
	*/
	//std::cerr << "BROKEN CHAIN COUNT: " << brokenChainCount << std::endl;
}

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_MSPLAZER_ALGORITHMS_H_
