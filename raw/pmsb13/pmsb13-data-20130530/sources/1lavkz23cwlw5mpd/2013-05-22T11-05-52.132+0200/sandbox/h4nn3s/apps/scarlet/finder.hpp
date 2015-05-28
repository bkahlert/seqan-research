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

#include "misc.hpp"

using namespace seqan;



//-----------------------------------------------------------------------------
//                  Finder
//-----------------------------------------------------------------------------

struct Match
{
    unsigned short qryId;
    unsigned short qryStart;
    unsigned short qryEnd;

    unsigned int   subjId; // many suffixes in subject-index
    unsigned short subjStart;
    unsigned short subjEnd;

    Match() :
        qryId(0), qryStart(0), qryEnd(0), subjId(0), subjStart(0), subjEnd(0)
    {
    }

    Match(Match const & m2)
    {
        qryId       = m2.qryId;
        qryStart    = m2.qryStart;
        qryEnd      = m2.qryEnd;
        subjId      = m2.subjId;
        subjStart   = m2.subjStart;
        subjEnd     = m2.subjEnd;
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
    String<Match>           matches;
    String<unsigned long>   seedRefs;  // mapping seed -> query
    String<unsigned int>    seedRanks; // mapping seed -> relative rank

    MatchManager()
    {
    }
};


template <typename TFinder>
inline void
onMatch(MatchManager & matchMan, TFinder const & finder)
{
    auto qryOccs = getOccurrences(back(finder.patternStack));
    auto subjOccs = getOccurrences(back(finder.textStack));

//     std::cout << length(qryOccs) << "\t" << length(subjOccs) << "\t"
//               << length(finder.patternStack) << std::endl;
//     for (auto const & s : finder.scoreStack)
//     {
//         for (auto const & ss : s)
//             std::cout << unsigned(ss) << ' ';
//         std::cout << std::endl;
//     }

    //TODO non-easy openmp
    for (auto const & qryOcc : qryOccs)
    {
        for (auto const & subjOcc : subjOccs)
        {
            Match m;
            m.qryId = getSeqNo(qryOcc);
            //TODO the following two can be removed, or not?
//             m.qryStart = getSeqOffset(qryOcc);
//             m.qryEnd = length(finder.patternStack) -2;
            m.subjId = getSeqNo(subjOcc);
            m.subjStart = getSeqOffset(subjOcc);
            m.subjEnd = getSeqOffset(subjOcc) + length(finder.textStack) -2;
            appendValue(matchMan.matches, m, Generous());
        }
    }
}

inline void
seedMatches2QueryMatches(MatchManager & mm, ScarletOptions const & options)
{
    // Update in place
    //TODO easy OpenMP
    for (unsigned long l = 0; l < length(mm.matches); ++l)
    {
        Match & m = value(mm.matches, l);
//         std::cout << "BEFORE:\n";
//         _printMatch(m);
        auto seedId = m.qryId;
        m.qryStart = mm.seedRanks[seedId] * options.seedLength;
        m.qryEnd = ((mm.seedRanks[seedId]+1) * options.seedLength) - 1;
        m.qryId = mm.seedRefs[seedId];
//         std::cout << "\nAFTER:\n";
//         _printMatch(m);
    }
}



inline void
joinAndFilterMatches(String<Match> & filteredMatches,
                     MatchManager & mm,
                     ScarletOptions const & options)
{
    unsigned char d = options.maxSeedDist;

    // sort
    std::cout << "Sorting matches..." << std::flush;
    double start = sysTime();
    std::sort(begin(mm.matches), end(mm.matches));
    double finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;

    // join and remove duplicates, filter too short matches
    std::cout << "Joining and filtering matches..." << std::flush;
    start = sysTime();

    String<unsigned char> seedsPerQuery;
    if (options.semiGlobal) // filter out matches that are too short
    {
        resize(seedsPerQuery,
               back(mm.seedRefs), // no of query sequences
               Exact());
        for (unsigned i = 0; i < length(mm.seedRefs); ++i)
            seedsPerQuery[mm.seedRefs[i]] = mm.seedRanks[i];
    }

    for (Match const & m : mm.matches)
    {
        if (length(filteredMatches) > 0)
        {
            // check if we can extend the last match with this one
            Match & fm = back(filteredMatches);
            if (// same sequence pair
                (m.qryId == fm.qryId) &&
                (m.subjId == fm.subjId) &&
                // and match beginnings overlap
                (((m.qryStart  >= _protectUnderflow(fm.qryStart, d)) &&
                  (m.qryStart  <= _protectOverflow( fm.qryEnd,   d)) &&
                  (m.subjStart >= _protectUnderflow(fm.subjStart,d)) &&
                  (m.subjStart <= _protectOverflow( fm.subjEnd,  d))
                 ) || // OR match endings overlap
                 ((m.qryEnd  >= _protectUnderflow(fm.qryStart, d)) &&
                  (m.qryEnd  <= _protectOverflow( fm.qryEnd,   d)) &&
                  (m.subjEnd >= _protectUnderflow(fm.subjStart,d)) &&
                  (m.subjEnd <= _protectOverflow( fm.subjEnd,  d))
               )))
            {
                fm.qryStart  = std::min(fm.qryStart, m.qryStart);
                fm.qryEnd    = std::max(fm.qryEnd,   m.qryEnd);
                fm.subjStart = std::min(fm.subjStart,m.subjStart);
                fm.subjEnd   = std::max(fm.subjEnd,  m.subjEnd);
                continue; // goto to next match
            }

            if (options.semiGlobal)
            {
                // if last match does not include all seeds
                if ((fm.qryEnd - fm.qryStart)
                      < (seedsPerQuery[fm.qryId] * options.seedLength))
                    // remove it
                    resize(filteredMatches, length(filteredMatches) - 1);
            }
        }
        // append the current match to the "output" list
        appendValue(filteredMatches, m, Generous());
    }
    finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "No of matches before joining and filtering: " << length(mm.matches)
              << ". After: " << length(filteredMatches) << "\n\n" << std::flush;

}

#endif // header guard
