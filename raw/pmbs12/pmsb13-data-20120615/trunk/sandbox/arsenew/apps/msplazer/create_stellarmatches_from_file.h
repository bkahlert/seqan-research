// ==========================================================================
//                      create_stellarmatches_from_file
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

#ifndef SANDBOX_KTRAPPE_APPS_MSPLAZER_CREATE_STELLARMATCHES_FROM_FILE_H_
#define SANDBOX_KTRAPPE_APPS_MSPLAZER_CREATE_STELLARMATCHES_FROM_FILE_H_

#include <iostream>
#include <fstream>

//#include <seqan/basic.h>
//#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.

//#include <seqan/stream.h>
#include <seqan/parse_lm.h>
#include "../../../../core/apps/stellar/stellar.h"
//#include "stellar_routines.h"

using namespace seqan;

//Creates a short Id out of a long one (i.e. it takes the prefix til the first white space)
template<typename TId>
void _getShortId(TId & longId, TId & shortId){

	clear(shortId);
	for (typename Position<TId>::Type i = 0; i < length(longId) && value(longId, i) > 32; ++i) {
    	    appendValue(shortId, value(longId, i));
    	}
}
//Creates a short Id out of a long one (i.e. it takes the prefix til the first white space)
template<typename TId>
void _getShortId(TId const & longId, TId & shortId){

	clear(shortId);
	for (typename Position<TId>::Type i = 0; i < length(longId) && value(longId, i) > 32; ++i) {
    	    appendValue(shortId, value(longId, i));
    	}
}

//Gets a Set of Ids and creates a set of short Ids from the longer ones
template<typename TId>
void _getShortId(StringSet<TId> & longIds, StringSet<TId> & shortIds){

	for (unsigned index = 0; index < length(longIds); ++index){
		TId & longId = longIds[index];
		TId sId;
		_getShortId(longId, sId);
		appendValue(shortIds, sId);
	}
}

//Takes the values from the localMatchStore and creates a Stellar match out of them
//The Stellar match is, sorted by read, appended to stQueryMatches
template < typename TSequence, typename TId> // typename TScoreAlloc > 
void _createStellarMatches(StringSet<TSequence> & queries, StringSet<TId> const & sQueryIds, StringSet<TSequence> & databases, StringSet<TId> const & sDBIds, StringSet<TId> const & databaseIds,
		StringSet<QueryMatches<StellarMatch<TSequence, TId> > > & stQueryMatches, LocalMatchStore<> & lmStore){//String <TScoreAlloc> & distScores, 

	typedef typename Infix<TSequence>::Type TInfix;
	typedef Segment<TInfix, InfixSegment> TSegment;
	typedef typename StellarMatch<TSequence, TId>::TAlign TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename LocalMatchStore<>::TPosition TPosition;


	for (unsigned i = 0; i < length(lmStore.matchStore); ++i) {

		TInfix dbInf;
	    TInfix queryInf;
		
		//Take match query ID and find the right index in queries
		//TPos iQuery = position of query/read sequence in queries
		unsigned iDB = -1;
		unsigned iQuery = -1;
		
		//Takes the short chromosome Id (from the Stellar match file) and looks up the corresponding long chromosome Id entry from the reference input file
		for(unsigned j = 0; j < length(sDBIds); ++j){
			if(lmStore.sequenceNameStore[lmStore.matchStore[i].subjectId] == sDBIds[j]){
				iDB = j;
				break;
			}
		}
		
		//Takes the short read Id (from the Stellar match file) and looks up the corresponding long read Id entry from the read input file
		for(unsigned j = 0; j < length(sQueryIds); ++j){
			if(lmStore.sequenceNameStore[lmStore.matchStore[i].queryId] == sQueryIds[j]){
				iQuery = j;
				break;
			}
		}
		//Sanity check for read and query Id: skippes entry if no corresponding entry in the input file could not be found, else creates StellarMatch object
		if(iDB == static_cast<unsigned>(-1) || iQuery == static_cast<unsigned>(-1)){
			std::cout << "Read or database does not exist for match: " << i << " subjectId: " << lmStore.sequenceNameStore[lmStore.matchStore[i].subjectId] << " queryId: " << lmStore.sequenceNameStore[lmStore.matchStore[i].queryId] << std::endl;
			std::cout << "Skipping entry" << std::endl;
		}
		else{

			//Checking orientation and swapping positions for reverse matches to apply them to stellar format
			bool orientation = true;
			if (lmStore.matchStore[i].subjectBeginPos > lmStore.matchStore[i].subjectEndPos){
				orientation = false;
				TPosition tmp = lmStore.matchStore[i].subjectBeginPos;
				lmStore.matchStore[i].subjectBeginPos = lmStore.matchStore[i].subjectEndPos;
				lmStore.matchStore[i].subjectEndPos = tmp;
			}
			//Computing infices for alignment rows for stellar matches
			queryInf = infix(queries[iQuery], lmStore.matchStore[i].queryBeginPos, lmStore.matchStore[i].queryEndPos);
			dbInf = infix(databases[iDB], lmStore.matchStore[i].subjectBeginPos, lmStore.matchStore[i].subjectEndPos);

			//Creating align object for stellar format
			TAlign localAlign;
			resize(rows(localAlign), 2);
    		setSource(row(localAlign, 0), host(dbInf));
	    	setSource(row(localAlign, 1), host(queryInf));
	       	TRow &row1 = row(localAlign,0);
			TRow &row2 = row(localAlign,1);


			// set begin and end positions of align
			setClippedBeginPosition(row2, lmStore.matchStore[i].queryBeginPos);
			setClippedBeginPosition(row1, lmStore.matchStore[i].subjectBeginPos);
			
			setBeginPosition(row1, 0);
			setBeginPosition(row2, 0);
			setClippedEndPosition(row1, lmStore.matchStore[i].subjectEndPos);
			setClippedEndPosition(row2, lmStore.matchStore[i].queryEndPos);
		
			unsigned gapIndex = 0;
			//Inserting gaps into rows according to cigar line
			if (length(lmStore.cigarStore) > lmStore.matchStore[i].id){
        	    String<CigarElement<> > const & cigar = lmStore.cigarStore[lmStore.matchStore[i].id];
            	for (unsigned j = 0; j < length(cigar); ++j){
                	//std::cout << cigar[j].count << cigar[j].operation;
					if(cigar[j].operation == 'I'){
						for(unsigned gap = 0; gap < cigar[j].count; ++gap){
							insertGap(row1,gapIndex);
							++gapIndex;
						}
					}
					else if(cigar[j].operation == 'D'){
						for(unsigned gap = 0; gap < cigar[j].count; ++gap){
							insertGap(row2,gapIndex);
							++gapIndex;
						}
					}
					else
						gapIndex += cigar[j].count;
				}
	        }
			//Create Stellar match and append it to stQueryMatches
			StellarMatch<TSequence, TId> match(localAlign, databaseIds[iDB], orientation);
			appendValue(stQueryMatches[iQuery].matches, match);
		}
		
	}
}

//Reads in a file with Stellar matches in gff format and creates StellarMatch object from the entries
template <typename TSequence, typename TId, typename TMatches > // typename TScoreAlloc >
void _getStellarMatchesFromFile(StringSet<TSequence> & queries, StringSet<TId> & queryIDs, StringSet<TSequence> & databases, StringSet<TId> & databaseIDs, 
		TMatches & stQueryMatches, CharString const & smFileName){
		//TMatches & stQueryMatches, String <TScoreAlloc> & distScores, CharString const & smFileName){

	//Allocator for short query Ids, needed bc Stellar only prints these short Ids to file
	StringSet<TId> sQueryIds;
	_getShortId(queryIDs, sQueryIds);
	//Allocator for short db Ids
	StringSet<TId> sDBIds;
	_getShortId(databaseIDs, sDBIds);
   
   	//Open file with Stellar matches	
	std::fstream inStreamMatches(toCString(smFileName), std::ios::in | std::ios::binary);
    if (!inStreamMatches.good()) {
    	std::cerr << "Could not open Stellar file " << smFileName << std::endl;
    }
	else{
		// Read local matches in GFF Stellar format.
    	RecordReader<std::fstream, SinglePass<> > recordReader(inStreamMatches);
	    LocalMatchStore<> lmStore;
	    unsigned i = 0;
    	while (!atEnd(recordReader))
	    {
        	int res = readRecord(lmStore, recordReader, StellarGff());
    	    if (res != 0)
            	std::cerr << "Invalid Stellar GFF record #" << i << '\n';
    	    i += 1;
	    }
		//Creating Stellar Matches from input
		resize(stQueryMatches, length(queries));

		//Infices/segments of db and query sequence, needed for stellar match
		typedef typename Infix<TSequence>::Type TInfix;
		typedef Segment<TInfix, InfixSegment> TSegment;

		_createStellarMatches(queries, sQueryIds, databases, sDBIds, databaseIDs, stQueryMatches, lmStore);
		//_createStellarMatches(queries, sQueryIds, databases, sDBIds, databaseIDs, stQueryMatches, distScores, lmStore);
	}
    //return 0;
}


#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_CREATE_STELLARMATCHES_FROM_FILE_H_
