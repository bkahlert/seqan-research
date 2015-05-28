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

// #include "options.hpp"
// #include "finder.hpp"


using namespace seqan;

// typedef StringSet<CharString, Owner<ConcatDirect<> > > TDbSeqs;
// typedef StringSet<CharString, Owner<ConcatDirect<> > > TQuerySeqs;
// 
// typedef Index<TDbSeqs,      IndexSa<> > TDbIndex;
// typedef Index<TQuerySeqs,   IndexSa<> > TQueryIndex;


struct BlastMatch_
{



};




// ============================================================================
// Forwards
// ============================================================================






// ============================================================================
// Metafunctions
// ============================================================================




// ============================================================================
// Functions for translation and retranslation
// ============================================================================



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

template <typename TString, typename TSpec, typename TFormat>
inline int
loadIds(StringSet<TString, TSpec > & ids,
              CharString const & path,
              TFormat const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    StringSet<CharString, TSpec > seqs;

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
// get plus-minus-range with bounds-checking for unsigned types
// ----------------------------------------------------------------------------

template <typename TNum, typename TNum2>
inline TNum
_protectUnderflow(const TNum n, const TNum2 s)
{
    const TNum r = n -s;
    return std::min(r, n);
}

template <typename TNum, typename TNum2>
inline TNum
_protectOverflow(const TNum n, const TNum2 s)
{
    const TNum r = n + s;
    return std::max(r, n);
}

template <typename TGaps>
inline bool
_startsWithGap(TGaps const & gaps)
{
    SEQAN_ASSERT(length(gaps._array) > 0);
    return (gaps._array[0] != 0);
}

template <typename TGaps>
inline int
_endsWithGap(TGaps const & gaps)
{
    SEQAN_ASSERT(length(gaps._array) > 0);
    if ((length(gaps._array)-1) % 2 == 1)
        return -1;
    return ((gaps._array[length(gaps._array)-1] != 0) ? 1 : 0);
}

template <typename TGaps, typename TSeq>
inline void
_prependNonGaps(TGaps & gaps, TSeq const & seq)
{
    if (_startsWithGap(gaps))
    {
        insertValue(gaps._array, 0, length(seq)); // new non-gap column
        insertValue(gaps._array, 0, 0); // empty gaps column
    }
    else
    {
        gaps._array[1] += length(seq);
    }

    insert(value(gaps._source), 0, seq);
    setBeginPosition(gaps, 0);
}

template <typename TGaps, typename TSeq>
inline void
_appendNonGaps(TGaps & gaps, TSeq const & seq)
{
    switch (_endsWithGap(gaps))
    {
        case -1:
        case  1:
            appendValue(gaps._array, length(seq)); // new non-gap column
            break;
        case 0:
            gaps._array[1] += length(seq);
    }
    append(value(gaps._source), seq);
    setEndPosition(gaps, length(value(gaps._source)));
}

#endif // header guard