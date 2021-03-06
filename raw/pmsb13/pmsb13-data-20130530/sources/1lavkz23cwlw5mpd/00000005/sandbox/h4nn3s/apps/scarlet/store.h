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
// store.h: contains types and definitions for storing sequences and indices
// ==========================================================================


// ============================================================================
// Forwards
// ============================================================================




// ============================================================================
// Tags, Classes, Enums
// ============================================================================


namespace seqan {
template <>
struct SAValue<TContigs>
{
    typedef Pair<unsigned int, unsigned int, Pack> Type;
};
}


typedef Index<TSubjects, IndexSa<> >                   TSubjectSA;

typedef MMap<>                                         TSubjectSaStringSpec;

namespace seqan {
template <>
struct Fibre<TGenomeSa, FibreSA>
{
    typedef String<SAValue<TContigs>::Type, TGenomeSaStringSpec>     Type;
};
}




// ============================================================================
// Metafunctions
// ============================================================================




// ============================================================================
// Functions
// ============================================================================