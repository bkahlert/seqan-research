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

#ifndef SANDBOX_KTRAPPE_APPS_MSPLAZER_MSPLAZER_H_
#define SANDBOX_KTRAPPE_APPS_MSPLAZER_MSPLAZER_H_

#define BREAKPOINT_DEBUG

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <seqan/misc/misc_cmdparser.h>

namespace seqan{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

///////////////////////////////////////////////////////////////////////////////
// MSplazer Options
struct MSplazerOptions {
	// i/o options
	CharString databaseFile;		// name of database file
	CharString queryFile;			// name of query file
	CharString outputFile;			// name of result file
	CharString samOutFile;		// name of sam result file
	CharString emptyReadsOutFile;	//	file with reads that have broken or empty chains
	CharString outDir;			// path for directory for graph files (.dot) and breakpoint files (.gff)
	CharString disabledQueriesFile;	// name of result file containing disabled queries
	CharString outputFormat;		// Possible formats: gff, text
	CharString jobName;			// job name, used for dot output
	CharString stellarInputFile;	// optional input file with stellar matches
	bool dotOut;

	//Transposition penalties and thresholds (criteria for inserting edges in read graphs)
	unsigned diffDBPen;				//Penalty for matches on different databases
	unsigned diffStrandPen;			//Pen. for matches on different strand of the same database
	unsigned diffOrderPen;			//Pen. for different order of matches within the same strand of one database with respect to the order within the read/query
	double simThresh;				//Allowed similarity between overlapping sequences
	int gapThresh;			//Allowed gap or distance between matches
	int initGapThresh;		//Maximal allowed start or end gap length
	unsigned breakpointPosRange;		//Allowed range of breakpoint positions
	unsigned support;
	

/*
	// main options for Stellar
	unsigned qGram;				// length of the q-grams
	double epsilon;				// maximal error rate
	int minLength;				// minimal length of an epsilon-match
	double xDrop;				// maximal x-drop
*/

	MSplazerOptions():
		databaseFile("reference.fa"),
		queryFile("reads.fa"),
		outputFile("msplazer.gff"),
		samOutFile("msplazer.sam"),
		emptyReadsOutFile("emptyReads.fa"),
		outDir(""),
		disabledQueriesFile("msplazer.disabled.fasta"),
		outputFormat("gff"),
		jobName(""),
		stellarInputFile(""),
		dotOut(false),
		diffDBPen(5),
		diffStrandPen(5),
		diffOrderPen(0),
		simThresh(0.5),
		gapThresh(10),
		initGapThresh(15),
		breakpointPosRange(6),
		support(2)
	{}

}; 

//////////////////////////////////////////////////////////////////////////////
// BreakPoint class: container for breakpoint information
template<typename TSequence_, typename TId_>
struct Breakpoint {
	typedef TSequence_							TSequence;
	typedef TId_								TId;
	typedef typename Position<TSequence>::Type	TPos;

	//Ids of the two sequences
	TId startSeqId;
	TId endSeqId;
	//Sequence orientation
	bool startSeqStrand;
	bool endSeqStrand;
	//Last position in start sequence and first position in end sequence
	TPos startSeqPos;
	TPos endSeqPos;
	TPos readStartPos;
	TPos readEndPos;
	//Counter of occurrences (read support)
	unsigned support;
	//Query Sequence Ids (queries/reads that support the breakpoint)
	StringSet<TId> supportIds;
	//SV type
	TId svtype;
	TSequence insertionSeq;
	bool revStrandDel;

	//For parameter testing: did bp result from an overlap (1) or a small insertion (0)
	bool overlapBP;


	Breakpoint (){}
	Breakpoint (TId const & sId, TId const & eId, bool const & sStrand, bool const & eStrand, TPos const & sPos, TPos const & ePos, TPos const & rsPos, TPos const & rePos):
		startSeqId(sId),
		endSeqId(eId),
		startSeqStrand(sStrand),
		endSeqStrand(eStrand),
		startSeqPos(sPos),
		endSeqPos(ePos),
		readStartPos(rsPos),
		readEndPos(rePos),
		support(1),
		svtype("svtype"),
		insertionSeq("NNNN"),
		revStrandDel(false)
	{}
	Breakpoint (TId const & sId, TId const & eId, bool const & sStrand, bool const & eStrand, TPos const & sPos, TPos const & ePos, TPos const & rsPos, TPos const & rePos, TId const & spId):
		startSeqId(sId),
		endSeqId(eId),
		startSeqStrand(sStrand),
		endSeqStrand(eStrand),
		startSeqPos(sPos),
		endSeqPos(ePos),
		readStartPos(rsPos),
		readEndPos(rePos),
		support(1),
		svtype("svtype"),
		insertionSeq("NNNN"),
		revStrandDel(false)
	{appendValue(supportIds, spId);}
	Breakpoint (TId const & sId, TId const & eId, bool const & sStrand, bool const & eStrand, TPos const & sPos, TPos const & ePos, TPos const & rsPos, TPos const & rePos, unsigned const & s, TId const & spId):
		startSeqId(sId),
		endSeqId(eId),
		startSeqStrand(sStrand),
   		endSeqStrand(eStrand),
		startSeqPos(sPos),
		endSeqPos(ePos),
		readStartPos(rsPos),
		readEndPos(rePos),
		support(s),
		svtype("svtype"),
		insertionSeq("NNNN"),
		revStrandDel(false)
	{appendValue(supportIds, spId);}
	Breakpoint (TId const & sId, TId const & eId, bool const & sStrand, bool const & eStrand, TPos const & sPos, TPos const & ePos, TPos const & rsPos, TPos const & rePos, unsigned const & s, StringSet<TId> const & spId):
		startSeqId(sId),
		endSeqId(eId),
   		startSeqStrand(sStrand),
		endSeqStrand(eStrand),
		startSeqPos(sPos),
		endSeqPos(ePos),
		readStartPos(rsPos),
		readEndPos(rePos),
		support(s),
		supportIds(spId),
		svtype("svtype"),
		insertionSeq("NNNN"),
		revStrandDel(false)
	{}

	//For testing
	Breakpoint (TId const & sId, TId const & eId, bool const & sStrand, bool const & eStrand, TPos const & sPos, TPos const & ePos, TPos const & rsPos, TPos const & rePos, TId const & spId, bool const & oBP):
		startSeqId(sId),
		endSeqId(eId),
		startSeqStrand(sStrand),
		endSeqStrand(eStrand),
		startSeqPos(sPos),
		endSeqPos(ePos),
		readStartPos(rsPos),
		readEndPos(rePos),
		support(1),
		svtype("svtype"),
		insertionSeq("NNNN"),
		revStrandDel(false),
		overlapBP(oBP)
	{appendValue(supportIds, spId);}

};

/////////////////////////////////////////////////////////////////////////////
//Sparse property map class: A property map that has only a few object or where most of the object would be empty
/**
.Class.SparsePropertyMap
..summary:Stores only a partial property map instead of one with many empty entries and a lookup table for the indices to the small property map.
..cat:Allocators
..signature:SparsePropertyMap<TValue, TPos>
..param.TValue:Type of stored objects.
..param.TPos:Type to store positions.
 
.Memfunc.SparsePropertyMap#SparsePropertyMap
..summary:Constructor
..signature:SparsePropertyMap<TValue,TPos> ()
..signature:SparsePropertyMap<TValue,TPos>  (TValueTable vt, TSlotLookupTable slt) 
..param.vt:Allocator for objects.
..param.slt:Allocator for positions.
..class:Class.SparsePropertyMap
.Memvar.SparsePropertyMap#valueTable
..summary:Allocator for objects.
..class:Class.SparsePropertyMap
.Memvar.SparsePropertyMap#slotLookupTable
..summary:Allocator for positions.
..class:Class.SparsePropertyMap
..include:msplazer.h
*/
template < typename TValue, typename TPos >
struct SparsePropertyMap {

	typedef String < TValue > 		TValueTable;
	typedef String < TPos >			TSlotLookupTable;

	TValueTable valueTable;
	TSlotLookupTable	slotLookupTable;

	SparsePropertyMap (){}
	//SparsePropertyMap (TValueTable vt, TSlotLookupTable slt) : valueTable(vt), slotLookupTable(slt) {}

};

template<typename TValue, typename TPos>
struct Value<SparsePropertyMap<TValue,TPos> > {
	typedef TValue Type;
};
template<typename TValue, typename TPos>
struct Value<SparsePropertyMap<TValue,TPos> const> {
	typedef TValue Type;
};

///////////////////////////////////////////////////////////////////////////////
// Container for storing chaining graph, matchDistanceScores, start and end vertex for one read
template < typename TGraph_, typename TVertexDescriptor_, typename TScoreAlloc_, typename TSparsePropertyMap_, typename TMatchAlloc_>
struct MSplazerChain {

	typedef TGraph_					TGraph;
	typedef TVertexDescriptor_ 		TVertexDescriptor;
	typedef TScoreAlloc_ 			TScoreAlloc;
	typedef TSparsePropertyMap_ 	TSparsePropertyMap;
	typedef typename Size<TGraph>::Type TGraphSize;
	typedef TMatchAlloc_		TMatchAlloc;

	//typedef TId_ 					TId;

    TGraph graph;						//Contains (Stellar)matches as vertices and edges between compatible matches
	TVertexDescriptor startVertex;		//Artificial start and end vertex (represents start and end of the read sequence)
	TVertexDescriptor endVertex;
	String<TVertexDescriptor> predMap;	//Predecessor and distance map for dagShortestPath
	String<TGraphSize> distMap;			//Distance map for shortest path
    TScoreAlloc matchDistanceScores;	//Distance scores of matches (edit distance of reference mapping)
	TSparsePropertyMap breakpoints;
	String<TMatchAlloc> bestChains;
	bool isEmpty;
	bool isPartial;
	//TId* readId;
	//TId* databaseId;

	//MSplazerChain (TGraph & g, TVertexDescriptor & sV, TVertexDescriptor & eV, TScoreAlloc & s) : graph(g),	startVertex(sV),
		//endVertex(eV), matchDistanceScores(s), isEmpty(false){}
	MSplazerChain (TScoreAlloc & _scores): matchDistanceScores(_scores), isEmpty(false), isPartial(false){}

};

struct Options
{
    bool showHelp;
    bool showVersion;
    int i;
    String<CharString> texts;
    
    Options()
    {
        // Set defaults.
        showHelp = true;
        i = 0;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/**
.Function.assignProperty:
..signature:assignProperty(spm, descr)
..param.spm: SparsePropertyMap.
...type:Class.SparsePropertyMap<TObject,TPos>
..param.descr: Vertex or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..include:msplazer.h
 */
template < typename TObject, typename TPos, typename TDescriptor >
inline void assignProperty(SparsePropertyMap<TObject, TPos> & spm, TDescriptor const & descr){
		assignProperty(spm.slotLookupTable, descr, -1);
}

/**
.Function.assignProperty:
..signature:assignProperty(spm, descr, val)
..param.spm: SparsePropertyMap.
...type:Class.SparsePropertyMap<TObject,TPos>
..param.descr: Vertex or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..param.val:The new value.
...remarks:Type of the new value must match the value type of the properties map object table.
.include:msplazer.h
 */
template < typename TObject, typename TPos, typename TDescriptor, typename TValue >
inline void assignProperty(SparsePropertyMap<TObject, TPos> & spm, TDescriptor const & descr, TValue const & val){
		assignProperty(spm.slotLookupTable, descr, length(spm.valueTable));
		appendValue(spm.valueTable, val);
}

/**
.Function.property:
..summary: Returns object at the given position. 
..signature:getProperty(spm, descr, obj)
..param.spm: SparsePropertyMap.
...type:Class.SparsePropertyMap<TObject,TPos>
..param.descr: Vertex or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..returns:Reference to the item in the property map.
..param.obj: Output parameter (object contained in sparse property map.
...remarks:Type must match the value type of the properties map object table. Assumes that there is an object at this position!
.include:msplazer.h
 */
template < typename TObject, typename TPos, typename TDescriptor >
inline typename Reference<TObject>::Type property(SparsePropertyMap<TObject, TPos> & spm, TDescriptor const & descr){

	TPos index = getProperty(spm.slotLookupTable, descr);
	//obj = getValue(spm.valueTable, index);
	return value(spm.valueTable, index);
}

template < typename TObject, typename TPos, typename TDescriptor >
inline typename Reference<TObject>::Type property(SparsePropertyMap<TObject, TPos> const & spm, TDescriptor const & descr){

	TPos index = getProperty(spm.slotLookupTable, descr);
	//obj = getValue(spm.valueTable, index);
	return value(spm.valueTable, index);
}


/**
.Function.getProperty:
..summary: Returns false if there is no object at the given position. Otherwise the object is written to the output parameter obj.
..signature:getProperty(spm, descr, obj)
..param.spm: SparsePropertyMap.
...type:Class.SparsePropertyMap<TObject,TPos>
..param.descr: Vertex or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..param.obj: Output parameter (object contained in sparse property map.
...remarks:Type must match the value type of the properties map object table.
.include:msplazer.h
 */
template < typename TObject, typename TPos, typename TDescriptor >
inline bool getProperty(SparsePropertyMap<TObject, TPos> & spm, TDescriptor const & descr, TObject & obj){

	TPos index = getProperty(spm.slotLookupTable, descr);
	if(index == static_cast<TPos>(-1))
		return 0;
	//obj = getValue(spm.valueTable, index);
	obj = value(spm.valueTable, index);
	return 1;
}

template < typename TObject, typename TPos, typename TDescriptor >
inline bool getProperty(SparsePropertyMap<TObject, TPos> const & spm, TDescriptor const & descr, TObject & obj){

	TPos index = getProperty(spm.slotLookupTable, descr);
	if(index == static_cast<TPos>(-1))
		return 0;
	//obj = getValue(spm.valueTable, index);
	obj = value(spm.valueTable, index);
	return 1;
}

/**
.Function.appendSupportId:
..signature:appendSupportId(bp, id)
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.id: Id to be appended.
...type:Class.TId
..include:msplazer.h
 */
template < typename TBreakpoint, typename TId >
inline void appendSupportId (TBreakpoint & bp, TId const & id)
{
	for(unsigned i = 0; i < length(bp.supportIds); ++i){
		if(bp.supportIds[i] == id){
			++bp.support;
			return;
		}
	}
	appendValue(bp.supportIds, id);
	++bp.support;
}

/**
.Function.appendSupportId:
..signature:appendSupportId(bp, id)
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.id: Id to be appended.
...type:TId
..include:msplazer.h
 */
template < typename TBreakpoint, typename TId >
inline void appendSupportId (TBreakpoint & bp, StringSet<TId> const & ids)
{

/*	
	typedef typename Iterator<StringSet<TId> >::Type TIterator;
	TIterator it = begin(ids);
	for(;!atEnd(it);goNext(it))
		appendSupportId(bp, *it);
*/		
	for(unsigned i = 0; i < length(ids); ++i)
		appendSupportId(bp, ids[i]);
}

/**
.Function.setSupport:
..signature:setSupport(bp, value)
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.value: New support value.
...type:TId
..include:msplazer.h
 */
template<typename TBreakpoint>
inline void setSupport (TBreakpoint & bp, unsigned const & value)
{
	assignValue(bp.support, value);
}

/**
.Function.setSVType:
..signature:setSVType(bp, type)
..summary:Sets the breakpoints svtype to "type".
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.type: TId.
...type:Class.Breakpoint.TId
..include:msplazer.h
 */
template<typename TBreakpoint, typename TId>
inline void setSVType(TBreakpoint & bp, TId const & type){
	bp.svtype = type;
}

/**
.Function.getSVType:
..signature:getSVType(bp)
..summary:Return the breakpoints svtype.
..param.bp: Breakpoint.
...type:Class.Breakpoint
..include:msplazer.h
 */
template<typename TSequence, typename TId>
inline TId getSVType(Breakpoint<TSequence, TId> & bp){
	return(bp.svtype);
}

/**
.Function.setSVType:
..signature:setSVType(bp)
..summary:Computes SV type of the beakpoint. Returns "true" for insertion and "false" otherwise.
..param.bp: Breakpoint.
...type:Class.Breakpoint
..include:msplazer.h
 */
template<typename TBreakpoint>
inline bool setSVType(TBreakpoint & bp){

	typedef typename TBreakpoint::TId TId;
	//clear(bp.svtype);
	//if insertion return 1; else return 0;
	if(bp.startSeqId != bp.endSeqId){
		bp.svtype = static_cast<TId>("translocation");
		//setSVType(bp, static_cast<TId>("translocation"));
		return 0;
	}
	if(bp.startSeqStrand != bp.endSeqStrand){
		bp.svtype = static_cast<TId>("inversion");
		//setSVType(bp, static_cast<TId>("inversion"));
		return 0;
	}
	if(bp.startSeqPos < bp.endSeqPos){
		bp.svtype = static_cast<TId>("deletion");
		//setSVType(bp, static_cast<TId>("deletion"));
		return 0;
	}
	if(bp.startSeqPos > bp.endSeqPos){
		std::swap(bp.startSeqPos, bp.endSeqPos);
		if(bp.startSeqStrand){
			bp.svtype = static_cast<TId>("translocation");
			return 0;
		}
		setSVType(bp, static_cast<TId>("deletion"));
		bp.revStrandDel = true;
		return 0;
	}
	bp.svtype = static_cast<TId>("insertion");
	//setSVType(bp, static_cast<TId>("insertion"));
	return 1;
}

/**
.Function.setInsertionSeq:
..signature:setInsertionSeq(bp, inSeq)
..summary:Computes SV type of the beakpoint. Returns true for insertion and false otherwise.
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.inSeq: Insertion sequence.
...type:TSequence
..include:msplazer.h
 */
template<typename TBreakpoint, typename TSequence>
inline void setInsertionSeq(TBreakpoint & bp, TSequence & inSeq)
{
	bp.insertionSeq = inSeq;
	//assignValue(bp.insertionSeq, inSeq);
}

/**
.Function.posInSameRange:
..signature:posInSameRange(pos1, pos2, range)
..param.pos1: First position to be compared.
...type:Class.Position
..param.pos2: Snd position to be compared.
...type:Class.Position
..param.range: Valid range for position difference.
...type:Class.Position
..include:msplazer.h
 */
template < typename TPos, typename TPosR >
inline bool posInSameRange (TPos const & pos1, TPos const & pos2, TPosR const & range)
{
	return abs(pos2 - pos1) < range;
}

/**
.Function.similarBreakpoints:
..summary:Tests two breakpoints for similarity, i.e. if they have the same sequence Ids and lie within a specified range.
..signature:similarBreakpoints(bp1, bp2i)
..param.bp1:First breakpoint to be compared.
...type:Class.Breakpoint
..param.bp2:Snd breakpoint to be compared.
...type:Class.Breakpoint
..include:seqan/msplazer.h
 */
template < typename TId, typename TPos >
inline bool similarBreakpoints	(Breakpoint<TId, TPos> const & bp1, Breakpoint<TId, TPos> const & bp2)
{
	bool sameSeqs = (bp1.startSeqId == bp2.startSeqId) && (bp1.endSeqId == bp2.endSeqId);
	bool samePosRange = posInSameRange(bp1.startSeqPos, bp2.startSeqPos) && posInSameRange(bp1.endSeqPos, bp2.endSeqPos);
	return sameSeqs && samePosRange;
}

/**
.Function.operator==:
..summary:Operator== implementation for Breakpoints. A Breakpoint has two Ids, two strands, two positions and a SV type, and can be compared according to them (in this order of priority).
..include:seqan/msplazer.h
 */
template < typename TId, typename TPos >
inline bool operator== (Breakpoint<TId, TPos> const & bp1, Breakpoint<TId, TPos> const & bp2)
{
	if(bp1.startSeqId != bp2.startSeqId)
		return false;   
	if(bp1.endSeqId != bp2.endSeqId)
		return false;
	/*
	if(bp1.startSeqStrand != bp2.startSeqStrand)
		return false;
	if(bp1.endSeqStrand != bp2.endSeqStrand)
		return false;
	*/	
	if(bp1.startSeqPos != bp2.startSeqPos)
		return false;
	if(bp1.endSeqPos != bp2.endSeqPos)
		return false;
	if(bp1.svtype != bp2.svtype)
		return false;
	if(bp1.svtype == "insertion" && bp2.svtype == "insertion"){
		//std::cerr << "Comparing two insertion breakpoints: " << length(bp1.insertionSeq) << " " << bp1.insertionSeq << " " << length(bp2.insertionSeq) << " " << bp2.insertionSeq << std::endl;
		//std::cerr << "returned: " << (length(bp1.insertionSeq) == length(bp2.insertionSeq)) << std::endl;
		return (length(bp1.insertionSeq) == length(bp2.insertionSeq));
	}
	return true;
}

/**
.Function.operator<:
..summary:Operator< implementation for Breakpoints. A Breakpoint has two Ids and two positions, and can be sorted according to them (in this order of priority).
..include:seqan/msplazer.h
 */
template < typename TId, typename TPos >
inline bool operator< (Breakpoint<TId, TPos> const & bp1, Breakpoint<TId, TPos> const & bp2)
{
	if(bp1.startSeqId != bp2.startSeqId)
		return (bp1.startSeqId < bp2.startSeqId);
	if(bp1.endSeqId != bp2.endSeqId)
		return (bp1.endSeqId < bp2.endSeqId);
	if(bp1.startSeqPos != bp2.startSeqPos)
		return (bp1.startSeqPos < bp2.startSeqPos);
	
	return (bp1.endSeqPos < bp2.endSeqPos);
}

/**
 * Ostream operator << for Breakpoint class
 */
template < typename TSequence, typename TId >
std::ostream& operator<< (std::ostream& out, Breakpoint< TSequence, TId > const & value)
{
	out << "Breakpoint: seq1 --> seq2; posInSeq1 --> posInSeq2; readPos1 --> readPos2 :" << std::endl;
	out << value.startSeqId << " ( " << value.startSeqStrand << " ) " << " --> " << value.endSeqId << " ( " << value.endSeqStrand << " ) " << std::endl;
	out << " ( " << value.startSeqPos + 1 << " ) --> ( " << value.endSeqPos + 1 << " ) " << std::endl;
	out << " ( " << value.readStartPos + 1 << " ) --> ( " << value.readEndPos + 1 << " ) " << std::endl;
	out << "SVType: " << value.svtype << " insertionSeq: " << value.insertionSeq << std::endl;
	out << "Support: " << value.support << " Ids: ";
	for(unsigned i = 0; i < length(value.supportIds); ++i)
		out << value.supportIds[i] << ", ";
	out << std::endl;
	return out;	
}

/**
 * Ostream operator << for StellarMatch
 */
template < typename TSequence, typename TId >
std::ostream& operator<< (std::ostream& out, StellarMatch<TSequence, TId> & match)
{
	out << "DB Id: " << match.id; 
		if(match.orientation)
			out << " + " << std::endl;
		else
			out << " - " << std::endl;
	out << "DB pos: " << match.begin1 << " ... " << match.end1 << std::endl;
	out << "Query pos: " << match.begin2 << " ... " << match.end2 << std::endl;

	if(!match.orientation)
		reverseComplement(infix(source(match.row1),match.begin1,match.end1));
	
	typedef typename StellarMatch<TSequence, TId>::TAlign TAlign;
	TAlign align;
	//resize(rows(align), 2);
	//assignSource makes a local copy of the source sequence (only the part of the row)
	//appendValue(row(align, 0), match.row1);
	//appendValue(row(align, 1), match.row2);
	appendValue(align.data_rows, match.row1);
	appendValue(align.data_rows, match.row2);
	out << align;
	out << std::endl;
	if(!match.orientation)
		reverseComplement(infix(source(match.row1),match.begin1,match.end1));
	return out;	
}
struct BreakpointFileHeader_;
typedef Tag<BreakpointFileHeader_> BreakpointFileHeader;

template <typename TSequence, typename TId>
inline void _streamPut(::std::FILE* target, Breakpoint<TSequence, TId> const & bp, unsigned const & counter){

	TId sId;
	//typedef typename Breakpoint<TSequence, TId>::TPos TPos;
	typedef typename Position<TSequence>::Type	TPos;
	_getShortId(bp.startSeqId, sId);

	_streamWrite(target, sId);
	_streamPut(target, '\t');
	_streamWrite(target, "MSplazer");
	_streamPut(target, '\t');
	_streamWrite(target, bp.svtype);
	_streamPut(target, '\t');
	_streamPutInt(target, bp.startSeqPos + 1);
	_streamPut(target, '\t');
	if(bp.svtype == "deletion" || bp.svtype == "inversion")
		_streamPutInt(target, bp.endSeqPos);
	else
		_streamPutInt(target, bp.startSeqPos + 1);
	_streamPut(target, '\t');
	_streamWrite(target, '.');
	_streamPut(target, '\t');
	if(bp.startSeqStrand)
		_streamWrite(target, '+');
	else
		_streamWrite(target, '-');
	_streamPut(target, '\t');
	_streamWrite(target, '.');
	_streamPut(target, '\t');

	_streamWrite(target, "ID=");	
	_streamPutInt(target, counter);
	_streamWrite(target, ';');
	if(bp.svtype == "insertion"){
		_streamWrite(target, "size=-");
		_streamPutInt(target, length(bp.insertionSeq));
		_streamWrite(target, ';');
		_streamWrite(target, "seq=");
		_streamWrite(target, bp.insertionSeq);
		_streamWrite(target, ';');
	}
	else if(bp.svtype == "deletion"){
		_streamWrite(target, "size=");
		_streamPutInt(target, static_cast<TPos>(bp.endSeqPos - bp.startSeqPos));
		_streamWrite(target, ';');
	}else{
		_streamWrite(target, "endChr=");
		_getShortId(bp.endSeqId, sId);
		_streamWrite(target, sId);
		_streamWrite(target, ';');
		//if(bp.svtype == "translocation"){
			_streamWrite(target, "endPos=");
			_streamPutInt(target, bp.endSeqPos);
			_streamWrite(target, ';');
		//}
		_streamWrite(target, "endStrand=");
		if(bp.endSeqStrand)
			_streamWrite(target, '+');
		else
			_streamWrite(target, '-');
		_streamWrite(target, ';');
	}
	_streamWrite(target, "support=");
	_streamPutInt(target, bp.support);
	_streamWrite(target, ';');
	_streamWrite(target, "supportIds=");
	for(unsigned i = 0; i < length(bp.supportIds); ++i){
		_streamWrite(target, bp.supportIds[i]);
		_streamWrite(target, ',');
	}
	_streamPut(target, '\n');
}
template <typename TSequence, typename TId>
inline void _streamWrite(::std::FILE* target, Breakpoint<TSequence, TId> const & bp, unsigned const & counter){
	_streamPut(target, bp, counter);
}
inline void _streamPut(::std::FILE* target, BreakpointFileHeader const &){

	_streamWrite(target, "StartSeqId");
	_streamPut(target, '\t');
	_streamWrite(target, "Label");
	_streamPut(target, '\t');
	_streamWrite(target, "SV type");
	_streamPut(target, '\t');
	_streamWrite(target, "sPos");
	_streamPut(target, '\t');
	_streamWrite(target, "sPos");
	_streamPut(target, '\t');
	_streamWrite(target, ".");
	_streamPut(target, '\t');
	_streamWrite(target, "+/-");
	_streamPut(target, '\t');
	_streamWrite(target, ".");
	_streamPut(target, '\t');
	_streamWrite(target, "Tags:ID=i;size=(-)indel;EndSeqId=;EndSeqPos=;EndSeqStrand=;Support=;SupportIds=;\n");
}
inline void _streamWrite(::std::FILE* target, BreakpointFileHeader const &){
	_streamPut(target, BreakpointFileHeader());
}

/**
.Function.BestChain:
..signature:write(msplChain, chain)
..param.msplChain: MSplazerChain object.
...type:Class.MSplazerChain
..param.chain: New chain.
...type:Class.MSplazerChain.TMatchAlloc
..include:seqan/msplazer.h
 */
template <typename TMSplazerChain, typename TMatchAlloc>
void insertBestChain(TMSplazerChain & mspChain, TMatchAlloc const & chain){

	appendValue(mspChain.bestChains, chain);
}

/**
.Function.insertBestChain:
..signature:write(msplChain, chain)
..param.msplChain: MSplazerChain object.
...type:Class.MSplazerChain
..param.chain: New chain.
...type:Class.MSplazerChain.TMatchAlloc
..include:seqan/msplazer.h
 */
template <typename TMSplazerChain, typename TMatchAlloc>
void insertBestChain(TMSplazerChain & mspChain, TMatchAlloc & chain){

	appendValue(mspChain.bestChains, chain);
}

//#include <sstream>

template <class T>
inline std::string toString (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

/**
.Tag.DotDrawingMSplazer
..cat:Input/Output
..summary:Switch to trigger drawing in dot format.
..value.DotDrawingMSplazer:Graphs in dot format.
..include:seqan/msplazer.h
*/

struct DotDrawingMSplazer_;
typedef Tag<DotDrawingMSplazer_> DotDrawingMSplazer;

/**
.Tag.DotDrawingMSplazerBestChain
..cat:Input/Output
..summary:Switch to trigger drawing in dot format.
..value.DotDrawingMSplazerBestChain:Best chain graphs in dot format.
..include:seqan/msplazer.h
*/

struct DotDrawingMSplazerBestChain_;
typedef Tag<DotDrawingMSplazer_> DotDrawingMSplazerBestchain;

/**
.Function.write:
..signature:write(file, msplazerchain, stellarmatches, tag)
..param.file:The file to write to.
...type:Class.Graph
..param.tag:A tag to select the output format.
...type:Tag.DotDrawingMSplazer
..include:seqan/msplazer.h
 */
template <typename TFile, typename TGraph, typename TVertexDescriptor, typename TScoreAlloc, typename TMatch, typename TBreakpointAlloc, typename TMatchAlloc > //, typename TId_ >
void 
write(TFile & file, 
	  MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, TBreakpointAlloc, TMatchAlloc > const & msplazerchain,
		TMatch const & queryMatches, unsigned const & queryLength,
	  	DotDrawingMSplazer const &) 
{

//IOREV _doc_ _batchreading_
	SEQAN_CHECKPOINT
	//typedef Graph<TSpec> TGraph;
	//typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Value<TBreakpointAlloc>::Type TBreakpoint;
	typedef typename TBreakpoint::TId TId;

	//_writeGraphType(file,g,DotDrawing());
	_streamWrite(file, "digraph G {\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Graph Attributes */\n");
	_streamWrite(file, "graph [rankdir = LR, clusterrank = local];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Node Attributes */\n");
	_streamWrite(file, "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Edge Attributes */\n");
	_streamWrite(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	_streamPut(file, '\n');

	_streamWrite(file, "/* Nodes */\n");
	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(msplazerchain.graph);
	unsigned i = 0;
	bool atEndV = false;
	for(;!atEnd(it);++it) {
	
		TId sId;	
		if(i < length(queryMatches)){
			_streamPutInt(file, *it);
			_streamWrite(file, " [label = \"");
			_streamWrite(file, "chr: ");
			_getShortId(queryMatches[i].id, sId);
			_streamWrite(file, sId);
			_streamWrite(file, "\\n db: ");
		//	_streamWrite(file, getProperty(nodeMap, *it));
		//std::cerr << " position: " << position(it) << std::endl;
			_streamPutInt(file, queryMatches[i].begin1 + 1);
			_streamWrite(file, "...");
			_streamPutInt(file, queryMatches[i].end1 + 1);
			_streamWrite(file, "  ");
			_streamWrite(file, (queryMatches[i].orientation ? '+' : '-'));
			_streamWrite(file, "\\n read: ");
			_streamPutInt(file, queryMatches[i].begin2 + 1);
			_streamWrite(file, "...");
			_streamPutInt(file, queryMatches[i].end2 + 1);
			_streamWrite(file, "\\n");
			_streamPutInt(file, i);

			_streamWrite(file, "\"];\n");
		}else if(!atEndV){
			_streamPutInt(file, *it);
			_streamWrite(file, " [label = \"start");

			_streamWrite(file, "\"];\n");
			atEndV = true;
		}else if(atEndV){
			_streamPutInt(file, *it);
			_streamWrite(file, " [label = \"end");
			_streamWrite(file, "\\n");
			_streamPutInt(file, queryLength + 1);

			_streamWrite(file, "\"];\n");

		}else
			std::cerr << "in writing dot: default vertex??" << std::endl;
		++i;
	}

	_streamPut(file, '\n');
	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(msplazerchain.graph);
	//bpC = 0;
	for(;!atEnd(itEd);++itEd) {
		TVertexDescriptor sc = sourceVertex(itEd);
		TVertexDescriptor tr = targetVertex(itEd);
		_streamPutInt(file, sc);
		_writeEdgeType(file, msplazerchain.graph, DotDrawing());
		_streamPutInt(file, tr);
		_streamWrite(file, " [label = \"");
		_streamPutInt(file, getCargo(*itEd));
		TBreakpoint bp;
		bool foundBP = getProperty(msplazerchain.breakpoints, value(itEd), bp);
		if(foundBP){
			_streamWrite(file, "*");
			//_streamPutInt(file, bpC);
			//++bpC;
		}
		_streamWrite(file, "\"];\n");
	}
	_streamPut(file, '\n');

	//_writeGraphFooter(file,g,DotDrawing()); //Doesnt do anything yet

	_streamWrite(file, "}\n");
}

/**TODO
.Function.write:
..signature:write(file, msplazerchain, stellarmatches, nodeMap, edgeMap, tag)
..param.graph:The graph to write out.
...type:Class.Graph
..param.nodeMap:A mapping from vertex descriptor to vertex label.
..param.edgeMap:A mapping from edge descriptor to edge label.
..param.tag:A tag to select the output format.
...type:Tag.SamMSplazer
..include:seqan/msplazer.h
 */

//////////////////////
//SAM Output
//TODO Use SamBamIO, see tutorial
template < typename TSequence, typename TId >
void _writeSamFile(CharString outfile, StringSet <QueryMatches< StellarMatch <TSequence, TId > > > & stellarMatches,  StringSet<TSequence > & databases, 
		StringSet<CharString> databaseIDs, StringSet<CharString> queryIDs){

	//CharString samOutFile2 = "msplazerStellarOut.sam";
	std::ofstream file2;
		file2.open(toCString(outfile), ::std::ios_base::out | ::std::ios_base::app);
        if (!file2.is_open()) {
                std::cerr << "Could not open output file " << outfile << std::endl;
        }else{

			typedef StellarMatch <TSequence, TId > TMatch;
			typedef typename Iterator<String<TMatch> >::Type TIterator;
    		typedef typename Row<typename TMatch::TAlign>::Type TMatchRow;
			typedef typename Size<TMatchRow>::Type TRowSize;
			TRowSize matchNumber, alignLen, alignDistance;

	      //  std::cerr << "length: " << length(stellarMatches[0].matches);


    	    file2    << "@HD\tVN:1.0\tSO:coordinate\n";
                for(unsigned i = 0; i < length(databaseIDs); ++i)
                        file2 << "@SQ\tSN:" << databaseIDs[i] << "\tLN:" << length(databases[i]) << "\n"; //matches[0].matches[0].id
                //TODO edit and flag CIGAR with clipped pos, mind strang direction
                for(unsigned i = 0; i < length(stellarMatches); ++i){

                        TIterator itStellarMatches = begin(stellarMatches[i].matches);
                        TIterator itEndStellarMatches = end(stellarMatches[i].matches);

        	        if (itStellarMatches != itEndStellarMatches) {
	                        stellarMatches[i].lengthAdjustment = _computeLengthAdjustment(length(source((*itStellarMatches).row1)), length(source((*itStellarMatches).row2)));
               		}

		      //  std::cerr << "test itStellarMatches: " << *itStellarMatches << " itEndStellarMatches: " << *itEndStellarMatches << std::endl;
                    while(itStellarMatches < itEndStellarMatches){

                        //for(unsigned stellarmatches = 0; stellarmatches < length(matches[i].matches); ++stellarmatches){
						if((*itStellarMatches).orientation){
	                                std::stringstream cigar, mutations;
        	                        _getCigarLine((*itStellarMatches).row1, (*itStellarMatches).row2, cigar, mutations);
  //            	                  std::cout << "cigar string: " << cigar.str() << std::endl;
//                      	          std::cout << "begin1: " << (*itStellarMatches).begin1 << std::endl;

							//Alignment Score using _analyzeAlignment(row0, row1, alignLen, matchNumber) from stellar_output.h
							//matchNumber and alignLen reset within _analyzeAlignment
							_analyzeAlignment((*itStellarMatches).row1, (*itStellarMatches).row2, alignLen, matchNumber);
							alignDistance = alignLen - matchNumber;
                                	file2 << queryIDs[i] << "\t"
	                                        << "0\t"
        	                                //<< databaseIDs[0] << "\t" //<< matches[i].matches[stellarmatches].id << "\t"
                	                        << (*itStellarMatches).id << "\t"
											//1-based position in DB/ref!
											<< "begin2: " << (*itStellarMatches).begin2 + 1 << "\t" //<< matches[i].matches[stellarmatches].begin2 << "\t"
											<< "end2: " << (*itStellarMatches).end2 + 1 << "\t" //<< matches[i].matches[stellarmatches].begin2 << "\t"
                                       		<< "begin1: " << (*itStellarMatches).begin1 + 1 << "\t" //<< matches[i].matches[stellarmatches].begin2 << "\t"
											<< "end1: " << (*itStellarMatches).end1 + 1 << "\t"
										<< ((*itStellarMatches).orientation ? '+' : '-') << "\t"
	                                        << "0\t"
        	                                << cigar.str() << "\t"
										<< "distance: " << alignDistance << "\t"
										<< "# matches: " << matchNumber << "\t"
										<< "alignLen: " << alignLen << "\t"
//                                      	  << "*\t0\t0\t"
	//                                        << "MatchSeq" << "\t"
        	                                << "*\n";
						}
                                ++itStellarMatches;//better goNext()?
					}

					reverseComplement(databases);
					itStellarMatches = begin(stellarMatches[i].matches);			

                    while(itStellarMatches < itEndStellarMatches){
                    //for(unsigned stellarmatches = 0; stellarmatches < length(matches[i].matches); ++stellarmatches){
						
						if(!(*itStellarMatches).orientation){

	                                std::stringstream cigar, mutations;
        	                        _getCigarLine((*itStellarMatches).row1, (*itStellarMatches).row2, cigar, mutations);
  //            	                  std::cout << "cigar string: " << cigar.str() << std::endl;
//                      	          std::cout << "begin1: " << (*itStellarMatches).begin1 << std::endl;
	
        	            //Alignment Score using _analyzeAlignment(row0, row1, alignLen, matchNumber) from stellar_output.h
                		//matchNumber and alignLen reset within _analyzeAlignment
                        _analyzeAlignment((*itStellarMatches).row1, (*itStellarMatches).row2, alignLen, matchNumber);
                        alignDistance = alignLen - matchNumber;
	                                file2 << queryIDs[i] << "\t"
        	                                << "0\t"
                	                        << databaseIDs[0] << "\t" //<< matches[i].matches[stellarmatches].id << "\t"
                        	                << (*itStellarMatches).id << "\t"
                                	        << "begin2: " << (*itStellarMatches).begin2 + 1 << "\t" //<< matches[i].matches[stellarmatches].begin2 << "\t"
											<< "end2: " << (*itStellarMatches).end2 + 1 << "\t" //<< matches[i].matches[stellarmatches].begin2 << "\t"
                                       		<< "begin1: " << (*itStellarMatches).begin1 + 1 << "\t" //<< matches[i].matches[stellarmatches].begin2 << "\t"
											<< "end1: " << (*itStellarMatches).end1 + 1 << "\t"
	                                    << ((*itStellarMatches).orientation ? '+' : '-') << "\t"
        	                                << "0\t"
                	                        << cigar.str() << "\t"
                        	            << "distance: " << alignDistance << "\t"
                                	    << "# matches: " << matchNumber << "\t"
                                        << "alignLen: " << alignLen << "\t"
	//                                        << "*\t0\t0\t"
//      	                                  << "MatchSeq" << "\t"
                	                        << "*\n";
						}
                                ++itStellarMatches;//better goNext()?
                     }
                     reverseComplement(databases);


                     file2 << "\n";
                }

        	file2.close();
		}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions from breakpoint.cpp (Anne-Katrin) until line 1395


// Edit distance match combination 
// places the breakpoint to the leftmost possible position if breakLeft==true
// to the rightmost position otherwise
template <typename TScore>
bool
findBestSplitPosition(String<Pair<TScore,int> > & maxColsL,
                                          String<Pair<TScore,int> > & maxColsR,
                                          int & rowPosL1,
                                          int rowPosL2,
                                          int & rowPosR1,
                                          int rowPosR2,
                                          int seq0Len,
                                          int & traceExt,
                                          bool breakLeft = true)
{

#ifdef BREAKPOINT_DEBUG
        //::std::cout << "findBestSplit\n";
#endif

        TScore maxSum = minValue<TScore>();
    int bestL = breakLeft ? rowPosL1 : rowPosL2;
    int bestR = breakLeft? rowPosR1 : rowPosR2;
        int bestTraceExtL = rowPosL1;
        int bestTraceExtR = rowPosL1;
        while (rowPosL1 <= rowPosL2 && rowPosR1 >= rowPosR2)
        {
                // this is to prevent same bases from being used in both prefix and suffix match 
                // this works, because we store the FIRST bestScore in each row
                if (!(maxColsL[rowPosL1].i2 + maxColsR[rowPosR1].i2 <= seq0Len))
                {
                        ++rowPosL1;
                        --rowPosR1;
                        continue;
                }
        // suboptimal
                if(maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1 > maxSum)
                {
                        maxSum = maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1;
                        bestL = rowPosL1;
                        bestR = rowPosR1;
                        bestTraceExtL = rowPosL1;
               }
                else if(maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1 == maxSum)
        {
            bestTraceExtR = rowPosL1;
            if (!breakLeft)
            {
                        maxSum = maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1;
                        bestL = rowPosL1;
                        bestR = rowPosR1;
            }
        }

                ++rowPosL1;
                --rowPosR1;
        }
        traceExt = bestTraceExtR - bestTraceExtL;
        rowPosL1 = bestL;
        rowPosR1 = bestR;

        return true;
}


// computes the score of a part of an alignment
template<typename TScoreValue,typename TScoreSpec,typename TAliSource,typename TAliSpec,typename TPosition, typename TSize>
TScoreValue
_getRefinedMatchScoreFromView(Score<TScoreValue,TScoreSpec> & score_type,
                 Align<TAliSource,TAliSpec> & segment,
                 TPosition view_pos,
                 TSize len)
{
SEQAN_CHECKPOINT
        typedef Align<TAliSource,TAliSpec> TAlign;
        typedef typename Row<TAlign>::Type TRow;
//      typedef typename Iterator<TRow,GapsIterator<ArrayGaps> >::Type TIterator;       
        typedef typename Iterator<TRow, Rooted>::Type TIterator;
        TIterator row0_it, row1_it;
        row0_it = iter(row(segment,0),view_pos);
        row1_it = iter(row(segment,1),view_pos);
        TSize i = 0;
        TScoreValue ret_score = 0;
    bool gap_is_open = false;
       while(i < len)
        {
                if(isGap(row1_it)||isGap(row0_it))
        {
            if(gap_is_open)     ret_score += scoreGapExtend(score_type);
            else
            {
                ret_score += scoreGapOpen(score_type);
                gap_is_open = true;
            }
        }
        else
        {
                        ret_score += score(score_type,getValue(row0_it),getValue(row1_it));
            gap_is_open = false;
        }
                ++i;
                ++row0_it;
                ++row1_it;
        }
        return ret_score;
}


// function combineMatchPair: looks for the optimal breakpoint position
// of two overlapping matches (match1 and match2)
// row1 specifies in which sequence the matches overlap (source sequence 0 or 1) in match1
// row2 specifies in which sequence the matches overlap (source sequence 0 or 1) in match2
// scoreType is the scoring scheme
// splitPos returns the best split position wrt the sequence specified with row1 and row2
// if trimMatches == true --> match1 and match2 are trimmed according to splitPos
// returns maxValue<TScoreValue>() if match combination doesn't work
// otherwise returns improvement in score induced by splitPos
template<typename TSequence, typename TScore,typename TSpec, typename TPosition>
typename Value<TScore>::Type
combineMatchPair(Align<TSequence,TSpec> &match1,
                 Align<TSequence,TSpec> &match2,
                 int row1,   // in which alignment row do the matches overlap? --> shortSeq
                 int row2,   // (i.e. read sequence for deletion, ref sequence for insertion)
                 bool reverse1,   
                 bool reverse2,   
                 TScore & scoreType,
                 TPosition &splitPos,
                 bool trimMatches = true) // directly trim matches according to best split position
{

    typedef typename Infix<TSequence>::Type TSequenceInf;
    typedef ModifiedString<TSequenceInf,ModReverse> TSequenceInfRev;
    typedef typename Value<TScore>::Type TScoreValue;

    // make sure row1 in match1 and row2 in match2 specify the same sequence
    SEQAN_ASSERT_EQ(getObjectId(source(row(match1,row1))),getObjectId(source(row(match2,row2))));

    // assume beginPos(match1) <= beginPos(match2)
    SEQAN_ASSERT_LEQ(clippedBeginPosition(row(match1,row1)),clippedBeginPosition(row(match2,row2)));

    // make sure match1 and match2 overlap
    if(clippedBeginPosition(row(match2,row2)) > clippedEndPosition(row(match1,row1)))
        return maxValue<TScoreValue>(); // --> match combination doesn't work

    // if one contains the other --> "discard" the smaller match
    if(clippedEndPosition(row(match2,row2)) <= clippedEndPosition(row(match1,row1)))
        return maxValue<TScoreValue>(); // --> match combination doesn't work

    //TPosition short1b = clippedBeginPosition(row(match1,row1));
    TPosition short1e = clippedEndPosition(row(match1,row1));
    TPosition short2b = clippedBeginPosition(row(match2,row2));
    //TPosition short2e = clippedEndPosition(row(match2,row2));

    TPosition long1b = clippedBeginPosition(row(match1,1-row1));
    TPosition long1e = clippedEndPosition(row(match1,1-row1));
 /*   TPosition long1eRev = long1e;
    if(reverse1)
    {
        TPosition long1btmp = long1b;
        long1b = length(source(row(match1,1-row1))) - long1e;
        long1e = length(source(row(match1,1-row1))) - long1btmp;
    }*/
    TPosition long2b = clippedBeginPosition(row(match2,1-row2));
    TPosition long2e = clippedEndPosition(row(match2,1-row2));
  /*  if(reverse2)
    {
        TPosition long2btmp = long2b;
        long2b = length(source(row(match2,1-row2))) - long2e;
        long2e = length(source(row(match2,1-row2))) - long2btmp;
    }*/

    // project short2b onto long sequence according to match1
    // --> begin position of first overlapping part on long seq
    TPosition short2bView1 = toViewPosition(row(match1,row1),short2b);
    TPosition longProj1 = toSourcePosition(row(match1,1-row1),short2bView1);
    //if(reverse1)
    //{
    //    longProj1 = length(source(row(match1,1-row1))) - longProj1;
    //}


    // project short1e onto long sequence according to match2
    // --> end position of second overlapping part on long seq
    TPosition short1eView2 = toViewPosition(row(match2,row2),short1e);
    TPosition longProj2 = toSourcePosition(row(match2,1-row2),short1eView2);
    //if(reverse2)
    //{
    //    longProj2 = length(source(row(match2,1-row2))) - longProj2;
    //}
    // first view position of match2

#ifdef BREAKPOINT_DEBUG
    //std::cout << "s1B=" << short1b << " s2B=" << short2b << " s1E=" << short1e << "s2E=" << short2e << std::endl;
    //std::cout << "l1B=" << long1b << " long1p=" << longProj1 << " s1e=" << long1e << "\t";
    //std::cout << "l2B=" << long2b << " long2p=" << longProj2 << " s2e=" << long2e << "\t";
#endif

    // overlapping part on short seq
    TSequenceInf shortOverlapSeq1 = infix(source(row(match1,row1)),short2b,short1e);
    TSequenceInfRev shortOverlapSeq2(infix(source(row(match1,row1)),short2b,short1e));

	//std::cerr << "shortOverlapSeq1: " << shortOverlapSeq1 << std::endl;
	//std::cerr << "shortOverlapSeq2: " << shortOverlapSeq2 << std::endl;

    //// process match1
    if(reverse1) reverseComplement(infix(source(row(match1,1-row1)),long1b,long1e));
    TSequenceInf longOverlapSeq1 = infix(source(row(match1,1-row1)),longProj1,long1e);


	//std::cerr << "longOverlapSeq1: " << longOverlapSeq1 << std::endl;

    //std::cout << match1;
    // count errors in overlapping part and set diagonals accordingly --> should not incorporate more errors
    // than are already in the alignments, then return error difference compared to before
    // last view position of match1
    TPosition viewEnd1 = _max(toViewPosition(row(match1,row1),short1e),toViewPosition(row(match1,1-row1),long1e/*Rev*/) );
    
    // get score of match1 in overlapping part
    int score1 = _getRefinedMatchScoreFromView(scoreType,match1,short2bView1,viewEnd1-short2bView1);

#ifdef BREAKPOINT_DEBUG
    //std::cout << "MatchScore of overlapping part in match1 = " << score1 << std::endl;
#endif
    // set diagonals for banded alignment
	//std::cerr << " " <<	-static_cast<int>(length(longOverlapSeq1)) << std::endl;
    int diag1L = _max(score1 - 1, -static_cast<int>(length(shortOverlapSeq1)));
    int diag2L = _min(-score1 + 1, static_cast<int>(length(longOverlapSeq1)));

    int minColNum = 0;
    // compute banded alignment of possible breakpoint region
    // rows in alignment matrix represent overlapSeq positions
    StringSet<TSequenceInf> strL;
    appendValue(strL,longOverlapSeq1);
    appendValue(strL,shortOverlapSeq1);
    String<Pair<TScoreValue,int> > maxColsL;
    Graph<Alignment<StringSet<TSequenceInf,Dependent<> >, void> > alignL(strL);

	/*
	std::cerr << "alignL: " << alignL << std::endl;
	for(unsigned i = 0; i < length(strL); ++i)
		std::cerr << "strL: " << strL[i] << std::endl;
	std::cerr << "diag1L " << diag1L << std::endl;
	std::cerr << "diag2L " << diag2L << std::endl;
	std::cerr << "maxColsL " << maxColsL << std::endl;
	std::cerr << "minColNum " << minColNum << std::endl;
	*/
    _globalAlignment(alignL,strL,scoreType,AlignConfig<false,false,false,false>(),diag1L,diag2L,maxColsL,minColNum,BandedNeedlemanWunsch());
    //std::cout << alignL;

    if(reverse1) reverseComplement(infix(source(row(match1,1-row1)),long1b,long1e));


    //// process match2
    if(reverse2) reverseComplement(infix(source(row(match2,1-row2)),long2b,long2e));
    TSequenceInfRev longOverlapSeq2(infix(source(row(match2,1-row2)),long2b,longProj2));

    //std::cout << match2;

    // first view pos of match2
    TPosition short2bView2 = 0;
    // get score of match2 in overlapping part
    int score2 = _getRefinedMatchScoreFromView(scoreType,match2,short2bView2,short1eView2);

#ifdef BREAKPOINT_DEBUG
    //std::cout << "MatchScore of overlapping part in match2 = " << score2 << std::endl;
#endif

    int diag1R = _max(score2 - 1, -(int)length(shortOverlapSeq2));
    int diag2R = _min(-score2 + 1, (int)length(longOverlapSeq2));

    StringSet<TSequenceInfRev> strR;
    appendValue(strR,longOverlapSeq2);
    appendValue(strR,shortOverlapSeq2);
    String<Pair<TScoreValue,int> > maxColsR;
    Graph<Alignment<StringSet<TSequenceInfRev,Dependent<> >, void> > alignR(strR);
    _globalAlignment(alignR,strR,scoreType,AlignConfig<false,false,false,false>(),diag1R,diag2R,maxColsR,minColNum,BandedNeedlemanWunsch());
    //std::cout << alignR;

    if(reverse2) reverseComplement(infix(source(row(match2,1-row2)),long2b,long2e));

    int rowPosL1 = 0;
        int rowPosR1 = length(shortOverlapSeq2);
        int rowPosL2 = length(shortOverlapSeq1);
        int rowPosR2 = 0;

    // compute best split position
    // seq0Len is only important in the case of a short deletion/splice to make sure that matches do not overlap on the other sequence (the 1-row1 and 1-row2 sequence = long seq)
    // only relevant of sequence is the same and order is match1 and match2 is also preserved on long seq
    int seq0Len = (getObjectId(source(row(match1,row1))) == getObjectId(source(row(match2,row2))) && longProj2 > longProj1) ? longProj2 - longProj1 : maxValue<int>();
        int traceExt = 0;
        findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,traceExt,true);

    splitPos = clippedBeginPosition(row(match2,row2)) + rowPosL1;

    if(trimMatches)
    {
        // TODO: trim matches according to splitPos
        // easy heuristic: trim according to view positions
        // exact: trim match1 according to short2b and append new trace (the one in alignL)
        //      + trim match2 according to short1e and pre-append new traces (from alignR)
        // for now: easy + faster version (the heuristic..)

        //// get relevant clipping positions
        TPosition projSplitPos1 = toSourcePosition(row(match1,1-row1),toViewPosition(row(match1,row1),splitPos));
//        int endGapsLong1 = projSplitPos1 > long1b ? toViewPosition(row(match2,1-row2),projSplitPos1) - toViewPosition(row(match2,1-row2),projSplitPos1-1) - 1 : 0;
        TPosition projSplitPos2 = toSourcePosition(row(match2,1-row2),toViewPosition(row(match2,row2),splitPos));
        //int numBeginGapsLong =  splitPos > short2b ? toViewPosition(row(match2,1-row2),projSplitPos2) - toViewPosition(row(match2,1-row2),projSplitPos2-1) - 1 : 0;
//        if(numBeginGapsLong > 0)
//            std::cout << "hier\n";
//        int numBeginGapsShort = splitPos > short2b ? toViewPosition(row(match2,row2),splitPos) - toViewPosition(row(match2,row2),splitPos-1) - 1 : 0;

        // TODO: trailing and leading gaps... not handled correctly yet!

        //// set positions on long sequence, by projecting splitPos
        setClippedEndPosition(row(match1,1-row1),projSplitPos1);
        setClippedBeginPosition(row(match2,1-row2),projSplitPos2);
        //setBeginPosition(row(match2,1-row2),numBeginGapsLong);  // the number of leading gaps

        // set positions on overlap sequencce
        setClippedEndPosition(row(match1,row1),splitPos);
        setClippedBeginPosition(row(match2,row2),splitPos);
        //setBeginPosition(row(match2,row2),0);    // gaps in the short sequence are consumed by the split gap
//       setBeginPosition(row(match2,row2),numBeginGapsShort);    // correct?

    }

    // return score difference:  overlapScoreNow - overlapScoreBefore

    return (maxColsR[rowPosR1].i1 + maxColsL[rowPosL1].i1) - (score1 + score2);

}

}//endof namespace
void
setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
    addVersionLine(parser, "0.1");
    
    addTitleLine(parser, "**********************");
    addTitleLine(parser, "* msplazer *");
    addTitleLine(parser, "**********************");
    addTitleLine(parser, "");
    addTitleLine(parser, "(c) 2011 by Your Name <your.email@example.net>");

    addUsageLine(parser, "msplazer [OPTIONS] TEXT+");
    
	addSection(parser, "Main Options");
	addOption(parser, CommandLineOption("i",  "integer",  "set an integer option", OptionType::Integer | OptionType::Label, options.i));
    
    requiredArguments(parser, 1);
}

inline void
_addVersion(CommandLineParser& parser) {
	::std::string rev = "$Revision: 8300 $";
	addVersionLine(parser, "Version 1.1 (February 5th 2011) SeqAn Revision: " + rev.substr(11, 4) + "");
}
///////////////////////////////////////////////////////////////////////////////
// Parses options from command line parser and writes them into options object
template<typename TParser, typename TOptions>
bool
_parseOptions(TParser & parser, TOptions & options, MSplazerOptions & msplazerOptions) {
//IOREV _notio_
    // i/o options
	getOptionValueShort(parser, 'd', options.databaseFile);
    getOptionValueShort(parser, 'q', options.queryFile);
	if (isSetShort(parser, 'i')) getOptionValueShort(parser, 'i', msplazerOptions.outDir);
	if (isSetShort(parser, 'j')) getOptionValueShort(parser, 'j', msplazerOptions.jobName);
    if (isSetShort(parser, 'o')) getOptionValueShort(parser, 'o', options.outputFile);
    if (isSetShort(parser, "od")) getOptionValueShort(parser, "od", options.disabledQueriesFile);
	if (isSetShort(parser, "of")) getOptionValueShort(parser, "of", options.outputFormat);
	if (isSetShort(parser, 'm')) getOptionValueShort(parser, 'm', msplazerOptions.stellarInputFile);
	if (isSetShort(parser, "do")) getOptionValueShort(parser, "do", msplazerOptions.dotOut);

	// main options
	
	if (isSetShort(parser, "tp")) getOptionValueShort(parser, "tp", msplazerOptions.diffDBPen);
	if (isSetShort(parser, "ip")) getOptionValueShort(parser, "ip", msplazerOptions.diffStrandPen);
	if (isSetShort(parser, "op")) getOptionValueShort(parser, "op", msplazerOptions.diffOrderPen);
	if (isSetShort(parser, "oth")) getOptionValueShort(parser, "oth", msplazerOptions.simThresh);
	if (isSetShort(parser, "gth")) getOptionValueShort(parser, "gth", msplazerOptions.gapThresh);
	if (isSetShort(parser, "ith")) getOptionValueShort(parser, "ith", msplazerOptions.initGapThresh);
	if (isSetShort(parser, "st")) getOptionValueShort(parser, "st", msplazerOptions.support);

	if (isSetLong(parser, "kmer")) getOptionValueLong(parser, "kmer", options.qGram);
    if (isSetLong(parser, "minLength")) getOptionValueLong(parser, "minLength", options.minLength);
	if (isSetShort(parser, 'e')) getOptionValueShort(parser, 'e', options.epsilon);
    if (isSetShort(parser, 'x')) getOptionValueShort(parser, 'x', options.xDrop);

	if (isSetShort(parser, 'f')) if (!isSetShort(parser, 'r')) options.reverse = false;
	if (isSetShort(parser, 'r')) if (!isSetShort(parser, 'f')) options.forward = false;

	if (isSetShort(parser, "vs")) getOptionValueShort(parser, "vs", options.fastOption);
	if (isSetShort(parser, "dt")) getOptionValueShort(parser, "dt", options.disableThresh);
	if (isSetShort(parser, 'n')) getOptionValueShort(parser, 'n', options.numMatches);
	if (isSetShort(parser, 's')) getOptionValueShort(parser, 's', options.compactThresh);
	if (isSetShort(parser, "rp")) getOptionValueShort(parser, "rp", options.maxRepeatPeriod);
	if (isSetShort(parser, "rl")) getOptionValueShort(parser, "rl", options.minRepeatLength);
	if (isSetShort(parser, 'a')) getOptionValueShort(parser, 'a', options.qgramAbundanceCut);

	if (isSetShort(parser, 'v')) options.verbose = 1;

	if (options.outputFormat != "gff" && options.outputFormat != "text") {
		std::cerr << "Invalid parameter value: Unknown output format." << std::endl;
		return 0;
	}

	if (options.fastOption != "exact" && options.fastOption != "bestLocal"
		 && options.fastOption != "bandedGlobal") {
		std::cerr << "Invalid parameter value: Unknown verification strategy." << std::endl;
		return 0;
	}

	if (isSetShort(parser, 'k') && options.qGram < 1) {
		std::cerr << "Invalid parameter value: Please choose a greater k-mer length." << std::endl;
		return 0;
	}

	if (isSetShort(parser, 'k') && options.qGram > 32) {
		std::cerr << "Invalid parameter value: Please choose a smaller k-mer length." << std::endl;
		return 0;
	}

	if (options.epsilon < 0.0) {
		std::cerr << "Invalid parameter value: Please choose a greater error rate." << std::endl;
		return 0;
	}

	if (options.epsilon > 0.25) {
		std::cerr << "Invalid parameter value: Please choose a smaller error rate." << std::endl;
		return 0;
	}
	if (isSetShort(parser, 'k') && options.qGram >= 1/options.epsilon) {
		std::cerr << "Invalid parameter value: Please choose q-gram length lower than 1/epsilon." << std::endl; 
		return 0;
	}

	if (options.qgramAbundanceCut > 1 || options.qgramAbundanceCut < 0) {
		std::cerr << "Invalid parameter value: Please choose a k-mer overabundance cut ration between 0 and 1.\n";
		return 0;
	}

	if (options.numMatches > options.compactThresh) {
		std::cerr << "Invalid parameter values: Please choose numMatches <= sortThresh." << std::endl;
		return 0;
	}
	return 1;
}

///////////////////////////////////////////////////////////////////////////////
// Set-Up of Command Line Parser
template<typename TParser>
void
_setParser(TParser & parser) {
	_addVersion(parser);

    addTitleLine(parser, "*******************************************");
	addTitleLine(parser, "*                MSplazer                 *");
	addTitleLine(parser, "* (c) Copyright 2011 by Kathrin Trappe    *");
	addTitleLine(parser, "*******************************************");

	addUsageLine(parser, "-d <FASTA sequence file> -q <FASTA sequence file> [Options]");

	addSection(parser, "MSplazer Options");	
	
	addSection(parser, "Non-optional Arguments:");
    addOption(parser, CommandLineOption('d', "database", "Fasta file containing the database sequences",
              (OptionType::String | OptionType::Mandatory)));
    addOption(parser, CommandLineOption('q', "query", "Fasta file containing the query sequences", 
              (OptionType::String | OptionType::Mandatory)));
	
	addSection(parser, "Main Options:");
	addOption(parser, CommandLineOption("tp", "transPen", "Interchromosomal translocation penalty", OptionType::Int, 5));
	addOption(parser, CommandLineOption("ip", "invPen", "Inversion penalty", OptionType::Int, 5));
	addOption(parser, CommandLineOption("op", "orderPen", "Intrachromosomal order change penalty", OptionType::Int, 0));
	addOption(parser, CommandLineOption("oth", "overlapThresh", "Allowed overlap between matches", OptionType::Double, "0.5"));
	addOption(parser, CommandLineOption("gth", "gapThresh", "Allowed gap length between matches", OptionType::Int, 10));
	addOption(parser, CommandLineOption("ith", "initGapThresh", "Allowed initial or ending gap length at begin and end of read", OptionType::Int, 15));
	addOption(parser, CommandLineOption("st", "support", "Number of supporting reads", OptionType::Int, 2));
	
	addSection(parser, "Input Options:");
	addOption(parser, CommandLineOption('m', "matchfile", "File of stellar matches", OptionType::String, ""));

	addSection(parser, "Output Options:");
	addOption(parser, CommandLineOption('i', "outdir", "Output directory", OptionType::String, ""));
	addOption(parser, CommandLineOption('j', "jobName", "Job/Queue name", OptionType::String, ""));
    addOption(parser, CommandLineOption("do", "dots", "Enable graph output in dot format", OptionType::Bool, "false"));
	
	addSection(parser, "Stellar Options");

	addSection(parser, "Main Options:");
    addOption(parser, CommandLineOption('e', "epsilon", "Maximal error rate (max 0.25)", OptionType::Double, "0.05"));
    addOption(parser, CommandLineOption('l', "minLength", "Minimal length of epsilon-matches", OptionType::Int, 100));
	addOption(parser, CommandLineOption('f', "forward", "Search only in forward strand of database",
		OptionType::Boolean, "both"));
	addOption(parser, CommandLineOption('r', "reverse", "Search only in reverse complement of database",
		OptionType::Boolean, "both"));
    addOption(parser, CommandLineOption('v', "verbose", "Verbosity mode.", OptionType::Bool, "false"));
    
	addSection(parser, "Filtering Options:");
    addOption(parser, CommandLineOption('k', "kmer", "Length of the q-grams (max 32)", OptionType::Int, 10));
    addOption(parser, CommandLineOption("rp", "repeatPeriod",
		"Maximal period of low complexity reapeats to be filtered", OptionType::Int, 1));
    addOption(parser, CommandLineOption("rl", "repeatLength",
		"Minimal length of low complexity reapeats to be filtered", OptionType::Int, 1000));
    addOption(parser, CommandLineOption('a', "abundanceCut",
		"k-mer overabundance cut ratio", OptionType::Double, "1"));

	addSection(parser, "Verification Options:");
    addOption(parser, CommandLineOption('x', "xDrop", "Maximal x-drop for extension", OptionType::Double, 5));
	addOption(parser, CommandLineOption("vs", "verification", "Verification strategy", OptionType::String, "exact"));
	addHelpLine(parser, "exact        = compute and extend all local alignments in SWIFT hits");
	addHelpLine(parser, "bestLocal    = compute and extend only best local alignment in SWIFT hits");
	addHelpLine(parser, "bandedGlobal = banded global alignment on SWIFT hits");
	addOption(parser, CommandLineOption("dt", "disableThresh",
		"Maximal number of verified matches before disabling verification", OptionType::Int));
	addHelpLine(parser, "for one query sequence (default infinity)");
	addOption(parser, CommandLineOption('n', "numMatches",
		"Maximal number of kept matches per query and database", OptionType::Int, 50));
	addHelpLine(parser, "If there are more matches, only the longest ones are kept.");
	addOption(parser, CommandLineOption('s', "sortThresh",
		"Number of matches triggering removal of duplicates", OptionType::Int, 500));
	addHelpLine(parser, "Choose a smaller value for saving space.");

	addSection(parser, "Output Options:");
    addOption(parser, CommandLineOption('o', "out", "Name of output file", OptionType::String, "stellar.gff"));
	addOption(parser, CommandLineOption("of", "outFormat", "Output format", OptionType::String, "gff"));
	addHelpLine(parser, "Possible formats: gff, text");
	addOption(parser, CommandLineOption("od", "outDisabled",
		"Name of output file containing disabled query sequences", OptionType::String));
	addHelpLine(parser, "(default stellar.disabled.fasta)");
}


int parseCommandLineAndCheck(Options & options,
                             CommandLineParser & parser,
                             int argc,
                             char const ** argv)
{
    bool stop = !parse(parser, argc, argv);
    if (stop)
        return 1;
    if (isSetLong(parser, "help")) {
        options.showHelp = true;
        return 0;
    }
    if (isSetLong(parser, "version")) {
        options.showVersion = true;
        return 0;
    }
    
    options.texts = getArgumentValues(parser);

	return 0;
}

int mainWithOptions(Options & options)
{
    typedef Iterator<String<CharString> >::Type TIterator;
    std::cout << "Non-option Arguments:" << std::endl;
    for (TIterator it = begin(options.texts); it != end(options.texts); ++it) {
        std::cout << "  " << *it << std::endl;
    }
    
    return 0;
}
#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_MSPLAZER_H_
