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

#ifndef SEQAN_SCARLET_STORE_H_
#define SEQAN_SCARLET_STORE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index.h>
#include <seqan/index_extras.h>

#include <seqan/blast.h>


using namespace seqan;

typedef StringSet<CharString, Owner<ConcatDirect<> > > TDbSeqs;
typedef StringSet<CharString, Owner<ConcatDirect<> > > TQuerySeqs;

typedef Index<TDbSeqs,      IndexSa<> > TDbIndex;
typedef Index<TQuerySeqs,   IndexSa<> > TQueryIndex;



// ============================================================================
// Forwards
// ============================================================================




// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Query
// ----------------------------------------------------------------------------

// template <BlastInputType>
// struct ScaSeq {};
// 
// template <typename TAlph, typename TSpec>
// struct ScaSeq<BlastInputType::Nucl>
// {
//     typedef 
//     const String<Dna5> seq;
// 
//     ScaSeq(String<TAlph, TSpec> const & in)
//         : seq(in)
//     {}
// };
// 
// template <typename TAlph, typename TSpec>
// struct ScaSeq<BlastInputType::Prot>
// {
//     const String<AminoAcid> seq;
// 
//     ScaSeq(String<TAlph, TSpec> const & in)
//         : seq(in)
//     {}
// };
// 
// template <typename TAlph, typename TSpec>
// struct ScaSeq<BlastInputType::TransNucl>
// {
//     const String<AminoAcid>   seq;
//     const unsigned long       origRank; // index of untranslated seq
//     const unsigned char       frameShift; // 0, 1 or 2
//     const bool                revComp; //reverseComplemented?
// 
//     ScaSeq(String<TAlph, TSpec> const & in,
//           const unsigned long inRank,
//           const unsigned char inShift,
//           const bool          inRevComp)
//         : seq(in), origRank(inRank), frameShift(inShift), revComp(inRevComp)
//     {}
// };

inline bool
isReverseComplemented(unsigned long l)
{
    if ((l % 6 == 0) || (l % 5 == 0) || (l % 4 == 0))
        return true;
}

inline unsigned char
getFrameShift(unsigned long l)
{
    return l % 3;
}

// ----------------------------------------------------------------------------
// Class ScaSeed
// ----------------------------------------------------------------------------


// std::tuple<unsigned long, unsigned char>


// template <typename TSpec>
// struct ScaSeed<ScaSeq<TSpec> >
// {
//     const unsigned long     parent; // index of the ScaSeq this belongs to
//     const unsigned char     n;     // this is n'th seed
//     const Infix<decltype(
// 
// };

// ============================================================================
// Metafunctions
// ============================================================================




// ============================================================================
// Functionshtop
// ============================================================================

template <typename TString, typename TSpec, typename TFormat>
inline int
loadSequences(StringSet<TString, TSpec > & seqs,
              CharString const & path,
              TFormat const & /*tag*/)
{
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




#endif // header guard