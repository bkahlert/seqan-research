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

struct MatchManager
{
    unsigned long numHits;

    MatchManager() :
        numHits(0)
    {
    }
} ;

//TODO -> EditDistance
typedef ext::Backtracking<HammingDistance, ext::BacktrackingSemiGlobal> BackSpec;
typedef ext::Finder<TDbIndex, TQueryIndex,
                    MatchManager, BackSpec> ScarletFinder;


template <typename TFinder>
void
onMatch(MatchManager & m, TFinder const & finder)
{
//     for (auto it = begin(finder.patternStack),
//          itEnd = end(finder.patternStack);
//         it != itEnd; ++it)
//         std::cout << position(it) << std::endl;
    ++m.numHits;
    std::cout << "pos(begin(VScoreStack)): "
              << position(begin(finder.scoreStack)) << '\t'
              << "pos(end(VScoreStack)): "
              << position(end(finder.scoreStack)) << '\n'
              << "pos(begin(patternStack)): "
              << position(begin(finder.patternStack)) << '\t'
              << "pos(end(patternStack)): "
              << position(end(finder.patternStack)) << '\n'
              << "pos(begin(textStack)): "
              << position(begin(finder.textStack)) << '\t'
              << "pos(end(textStack)): "
              << position(end(finder.textStack)) << "\n\n";

//     std::cout << "Found something\n";
}






#endif // header guard