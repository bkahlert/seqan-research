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

#ifndef SEQAN_SCARLET_MISC_H_
#define SEQAN_SCARLET_MISC_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index.h>
#include <seqan/index_extras.h>

#include <seqan/blast.h>

#include "options.hpp"


using namespace seqan;

// typedef StringSet<CharString, Owner<ConcatDirect<> > > TDbSeqs;
// typedef StringSet<CharString, Owner<ConcatDirect<> > > TQuerySeqs;
// 
// typedef Index<TDbSeqs,      IndexSa<> > TDbIndex;
// typedef Index<TQuerySeqs,   IndexSa<> > TQueryIndex;



// ============================================================================
// Forwards
// ============================================================================






// ============================================================================
// Metafunctions
// ============================================================================




// ============================================================================
// Functions for translation and retranslation
// ============================================================================

template <typename TNum>
inline TNum
getTruePos(TNum pos, bool isDb, ScarletOptions const & options)
{
    if (isDb) // is subject sequence
        switch (options.blastProg)
        {
            case BlastFormatOptions::TBlastN :
            case BlastFormatOptions::TBlastX : 
                return pos / 6; // 6 frame translated nucleotide
            default:
                break;
        }

    // is query sequence
    switch (options.blastProg)
    {
        case BlastFormatOptions::BlastN :
            if (options.revComp)
                return pos / 2; // has reverse complements
            else
                break;
        case BlastFormatOptions::BlastX :
        case BlastFormatOptions::TBlastX :
            return pos / 6; // 6 frame translated nucleotide
        default:
            break;
    }
    return pos;
}

template <typename TNum>
inline bool
isReverseComplemented(TNum pos, bool isDb, ScarletOptions const & options)
{
    if (isDb) // is subject sequence
        switch (options.blastProg)
        {
            case BlastFormatOptions::TBlastN :
            case BlastFormatOptions::TBlastX :
                // 6 frame translated nucleotide
                return ((pos % 6 == 0) || (pos % 5 == 0) || (pos % 4 == 0));
            default:
                break;
        }

    // is query sequence
    switch (options.blastProg)
    {
        case BlastFormatOptions::BlastN :
            if (options.revComp)
                return (pos % 2 == 0); // has reverse complements
            else
                break;
        case BlastFormatOptions::BlastX :
        case BlastFormatOptions::TBlastX :
            // 6 frame translated nucleotide
            return ((pos % 6 == 0) || (pos % 5 == 0) || (pos % 4 == 0));
        default:
            break;
    }
    return false;
}

template <typename TNum>
inline unsigned char
getFrameShift(TNum pos, bool isDb, ScarletOptions const & options)
{
    if (isDb) // is subject sequence
        switch (options.blastProg)
        {
            case BlastFormatOptions::TBlastN :
            case BlastFormatOptions::TBlastX :
                // 6 frame translated nucleotide
                return (pos % 3);
            default:
                break;
        }

    // is query sequence
    switch (options.blastProg)
    {
        case BlastFormatOptions::BlastX :
        case BlastFormatOptions::TBlastX :
            // 6 frame translated nucleotide
            return (pos % 3);
        default:
            break;
    }
    return 0;
}


// ----------------------------------------------------------------------------
// Generic Sequence loading
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TFormat>
inline int
loadSequences(StringSet<TString, TSpec > & seqs,
              CharString const & path,
              TFormat const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    StringSet<CharString, TSpec > ids;

    int res = read2(ids, seqs, reader, TFormat());
    if (res)
        std::cerr << "Error : " << res << "\n";

    stream.close();
    return res;
}


// ----------------------------------------------------------------------------
// truncate sequences
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline void
_debug_shorten(StringSet<TString, TSpec > & seqs, unsigned const len)
{
    StringSet<TString, TSpec > copySeqs;
    reserve(copySeqs.concat, length(seqs)*len, Exact());

    for (TString const & s : seqs)
        if (length(s) >= len)
            appendValue(copySeqs, prefix(s, len), Exact());

    clear(seqs);
    reserve(seqs.concat, length(copySeqs)*len, Exact());
    for (TString const & s : copySeqs)
        appendValue(seqs, s);
}

// ----------------------------------------------------------------------------
// get plus-minus-range with bounds-checking
// ----------------------------------------------------------------------------

template <typename TNum, typename TNum2>
inline TNum
_protectUnderflow(const TNum n, const TNum2 s)
{
    return (n - s < n) ?  n - s : n;
}

template <typename TNum, typename TNum2>
inline TNum
_protectOverflow(const TNum n, const TNum2 s)
{
    return (n + s > n) ?  n + s : n;
}



#endif // header guard