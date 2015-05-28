// ==========================================================================
//                                  scarlet
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
// All rights reserved.
//
// This file is part of Scarlet.
//
// Scarlet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Scarlet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Scarlet.  If not, see <http://www.gnu.org/licenses/>.*/
// ==========================================================================
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// finder.h: Main File for the indexer application
// ==========================================================================

#ifndef SEQAN_SCARLET_FINDER_H_
#define SEQAN_SCARLET_FINDER_H_

#include "store.hpp"

using namespace seqan;


//-----------------------------------------------------------------------------
//                  Finder
//-----------------------------------------------------------------------------

struct Match
{
    unsigned short qryId;
    unsigned short qryStart;  // = 0
    unsigned short qryEnd;    // = end

    unsigned short subjId;
    unsigned short subjStart;
    unsigned short subjEnd;

    Match() :
        qryId(0), qryStart(0), qryEnd(0), subjId(0), subjStart(0), subjEnd(0)
    {
    }

    inline bool operator< (Match const & m2) const
    {
        if (qryId != m2.qryId)
            return qryId < m2.qryId;
        if (subjId != m2.subjId)
            return subjId < m2.subjId;
        if (qryStart != m2.qryStart)
            return qryStart < m2.qryStart;
        if (subjStart != m2.subjStart)
            return subjStart < m2.subjStart;
        if (qryEnd != m2.qryEnd)
            return qryEnd < m2.qryEnd;
        if (subjEnd != m2.subjEnd)
            return subjEnd < m2.subjEnd;
        return false;
    }
};

inline void
_printMatch(Match const & m)
{
    std::cout << "MATCH  Query " << m.qryId
              << "(" << m.qryStart << ", " << m.qryEnd
              << ")   on Subject "<< m.subjId
              << "(" << m.subjStart << ", " << m.subjEnd
              << ")" <<  std::endl << std::flush;
}

struct MatchManager
{
    String<Match> matches;

    MatchManager()
    {
    }
} ;

//TODO -> EditDistance
typedef ext::Backtracking<EditDistance, ext::BacktrackingSemiGlobal> BackSpec;
typedef ext::Finder<TDbIndex, TQueryIndex,
                    MatchManager, BackSpec> ScarletFinder;


template <typename TFinder>
inline void
onMatch(MatchManager & matchMan, TFinder const & finder)
{
    auto qryOccs = getOccurrences(back(finder.patternStack));
    auto subjOccs = getOccurrences(back(finder.textStack));

    std::cout << length(qryOccs) << "\t" << length(subjOccs) << "\t"
              << length(finder.patternStack) << std::endl;
    for (auto const & s : finder.scoreStack)
    {
        for (auto const & ss : s)
            std::cout << unsigned(ss) << ' ';
        std::cout << std::endl;
    }

    for (auto const & qryOcc : qryOccs)
    {
        for (auto const & subjOcc : subjOccs)
        {
            Match m;
            m.qryId = getSeqNo(qryOcc);
            m.qryStart = getSeqOffset(qryOcc);
            m.qryEnd = length(finder.patternStack) -2;
            m.subjId = getSeqNo(subjOcc);
            m.subjStart = getSeqOffset(subjOcc);
            m.subjEnd = length(finder.textStack) -2;
            appendValue(matchMan.matches, m, Generous());
        }
    }
}






#endif // header guard
