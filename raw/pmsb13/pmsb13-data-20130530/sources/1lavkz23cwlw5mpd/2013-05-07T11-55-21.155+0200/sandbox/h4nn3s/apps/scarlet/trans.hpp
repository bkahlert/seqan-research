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
// trans.hpp: Functions and Metafunctions for translation
// ==========================================================================

#ifndef SEQAN_SCARLET_TRANS_H_
#define SEQAN_SCARLET_TRANS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index.h>
#include <seqan/index_extras.h>

#include <seqan/blast.h>

#include "options.hpp"
#include "finder.hpp"

using namespace seqan;


char translationTable[4][4][4] =
{
    { // a??
        { 'K', 'N', 'K', 'N' }, // aa?
        { 'T', 'T', 'T', 'T' }, // ac?
        { 'R', 'S', 'R', 'S' }, // ag?
        { 'I', 'I', 'M', 'I' }  // au?
    },
    { // c??
        { 'Q', 'H', 'Q', 'H' }, // ca?
        { 'P', 'P', 'P', 'P' }, // cc?
        { 'R', 'R', 'R', 'R' }, // cg?
        { 'L', 'L', 'L', 'L' }  // cu?
    },
    { // g??
        { 'E', 'D', 'E', 'D' }, // ga?
        { 'A', 'A', 'A', 'A' }, // gc?
        { 'G', 'G', 'G', 'G' }, // gg?
        { 'V', 'V', 'V', 'V' }  // gu?
    },
    { // u??
        { '*', 'Y', '*', 'Y' }, // ua?
        { 'S', 'S', 'S', 'S' }, // uc?
        { '*', 'C', 'W', 'C' }, // ug?
        { 'L', 'F', 'L', 'F' }  // uu?
    }
};


// ----------------------------------------------------------------------------
// append6FrameTranslation()
// ----------------------------------------------------------------------------

inline void
appendSixFrameTranslation(StringSet<AminoAcid, Owner<ConcatDirect<>>> & target,
                          String<Dna5> & source) // not const and will be modified!
{
    for (unsigned char c : { 0, 1})
    {
        for (unsigned char cc : { 0, 1, 2 } )
        {
            appendValue(target, String<AminoAcid>());
//             reserve(back(target), (length(source) - cc) / 3, Exact());
            for (unsigned int i = cc; i < length(source) - 3; i+=3)
            {
                if ((ordValue(value(source, i  ) == 4) ||
                    (ordValue(value(source, i+1) == 4) || // 'N' encountered
                    (ordValue(value(source, i+2) == 4) )
                    appendValue(back(target), 'X', Exact());
                else
                    appendValue(back(target),
                                translationTable[ordValue(value(source, i  ))]
                                                [ordValue(value(source, i+1))]
                                                [ordValue(value(source, i+2))],
                                Exact());
            }
        }
        reverseComplement(source);
    }
}
// ----------------------------------------------------------------------------
// getTrueQryId()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryId(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.qryId)
{
    return m.qryId;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryId(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastN,g> const & /*tag*/)
-> decltype(m.qryId)
{
    return (options.revComp) ? (m.qryId / 2) : m.qryId;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryId(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastX,g> const & /*tag*/)
-> decltype(m.qryId)
{
    return m.qryId / 3;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryId(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.qryId)
{
    return m.qryId / 3;
}

// ----------------------------------------------------------------------------
// getTrueSubjId()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjId(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.subjId)
{
    return m.subjId;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjId(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastN,g> const & /*tag*/)
-> decltype(m.subjId)
{
    return m.subjId / 3;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjId(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.subjId)
{
    return m.subjId / 3;
}

// ----------------------------------------------------------------------------
// qryIsReverseComplemented()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr bool
qryIsReverseComplemented(Match const & m,
                      ScarletOptions const & options,
                      BlastFormat<mf,p,g> const & /*tag*/)
{
    return false;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
qryIsReverseComplemented(Match const & m,
                      ScarletOptions const & options,
                      BlastFormat<mf,
                                  BlastFormatOptions::BlastN,
                                  g> const & /*tag*/)
{
    return (options.revComp) ? (m.qryId % 2 == 0) : false;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
qryIsReverseComplemented(Match const & m,
                      ScarletOptions const & options,
                      BlastFormat<mf,
                                  BlastFormatOptions::BlastX,
                                  g> const & /*tag*/)
{
    return ((m.qryId+1 % 4 == 0) ||
            (m.qryId+1 % 5 == 0) ||
            (m.qryId+1 % 6 == 0));
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
qryIsReverseComplemented(Match const & m,
                      ScarletOptions const & options,
                      BlastFormat<mf,
                                  BlastFormatOptions::TBlastX,
                                  g> const & /*tag*/)
{
    return ((m.qryId+1 % 4 == 0) ||
            (m.qryId+1 % 5 == 0) ||
            (m.qryId+1 % 6 == 0));
}

// ----------------------------------------------------------------------------
// subjIsReverseComplemented()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr bool
subjIsReverseComplemented(Match const & m,
                      ScarletOptions const & options,
                      BlastFormat<mf,p,g> const & /*tag*/)
{
    return false;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
subjIsReverseComplemented(Match const & m,
                      ScarletOptions const & options,
                      BlastFormat<mf,
                                  BlastFormatOptions::TBlastN,
                                  g> const & /*tag*/)
{
    return ((m.subjId+1 % 4 == 0) ||
            (m.subjId+1 % 5 == 0) ||
            (m.subjId+1 % 6 == 0));
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
subjIsReverseComplemented(Match const & m,
                      ScarletOptions const & options,
                      BlastFormat<mf,
                                  BlastFormatOptions::TBlastX,
                                  g> const & /*tag*/)
{
    return ((m.subjId+1 % 4 == 0) ||
            (m.subjId+1 % 5 == 0) ||
            (m.subjId+1 % 6 == 0));
}

// ----------------------------------------------------------------------------
// getQryFrameShift()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getQryFrameShift(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,p,g> const & /*tag*/)
{
    return 0;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getQryFrameShift(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastX,g> const & /*tag*/)
{
    return m.qryId % 3;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getQryFrameShift(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
{
    return m.qryId % 3;
}

// ----------------------------------------------------------------------------
// getSubjFrameShift()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getSubjFrameShift(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,p,g> const & /*tag*/)
{
    return 0;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getSubjFrameShift(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastN,g> const & /*tag*/)
{
    return m.subjId % 3;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr unsigned char
getSubjFrameShift(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
{
    return m.subjId % 3;
}


// ----------------------------------------------------------------------------
// getTrueQryStartPos()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryStartPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.qryStart)
{
    return m.qryStart;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryStartPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastN,g> const & /*tag*/)
-> decltype(m.qryStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::BlastN,g> TFormat;
    return (qryIsReverseComplemented(m, options, TFormat()))
                ? m.qryEnd
                : m.qryStart;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryStartPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastX,g> const & /*tag*/)
-> decltype(m.qryStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::BlastX,g> TFormat;
    return (qryIsReverseComplemented(m, options, TFormat()))
            ? m.qryEnd   * 3 - getQryFrameShift(m, options, TFormat())
            : m.qryStart * 3 + getQryFrameShift(m, options, TFormat());
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryStartPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.qryStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastX,g> TFormat;
    return (qryIsReverseComplemented(m, options, TFormat()))
            ? (m.qryEnd   * 3 - getQryFrameShift(m, options, TFormat()))
            : (m.qryStart * 3 + getQryFrameShift(m, options, TFormat()));
}


// ----------------------------------------------------------------------------
// getTrueQryEndPos()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryEndPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.qryEnd)
{
    return m.qryEnd;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryEndPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastN,g> const & /*tag*/)
-> decltype(m.qryEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::BlastN,g> TFormat;
    return (isReverseComplemented(mf, options, TFormat()))
                ? m.qryStart
                : m.qryEnd;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryEndPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::BlastX,g> const & /*tag*/)
-> decltype(m.qryEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::BlastX,g> TFormat;
    return (!qryIsReverseComplemented(m, options, TFormat()))
            ? (m.qryEnd   * 3 - getQryFrameShift(m, options, TFormat()))
            : (m.qryStart * 3 + getQryFrameShift(m, options, TFormat()));
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueQryEndPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.qryEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastX,g> TFormat;
    return (!qryIsReverseComplemented(m, options, TFormat()))
            ? (m.qryEnd   * 3 - getQryFrameShift(m, options, TFormat()))
            : (m.qryStart * 3 + getQryFrameShift(m, options, TFormat()));
}

// ----------------------------------------------------------------------------
// getTrueSubjStartPos()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjStartPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.subjStart)
{
    return m.subjStart;
}


template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjStartPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastN,g> const & /*tag*/)
-> decltype(m.subjStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastN,g> TFormat;
    return (subjIsReverseComplemented(m, options, TFormat()))
            ? (m.subjEnd   * 3 - getQryFrameShift(m, options, TFormat()))
            : (m.subjStart * 3 + getQryFrameShift(m, options, TFormat()));
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjStartPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.subjStart)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastX,g> TFormat;
    return (subjIsReverseComplemented(m, options, TFormat()))
            ? (m.subjEnd   * 3 - getQryFrameShift(m, options, TFormat()))
            : (m.subjStart * 3 + getQryFrameShift(m, options, TFormat()));
}


// ----------------------------------------------------------------------------
// getTrueSubjEndPos()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjEndPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,p,g> const & /*tag*/)
-> decltype(m.subjEnd)
{
    return m.subjEnd;
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjEndPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastN,g> const & /*tag*/)
-> decltype(m.subjEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastN,g> TFormat;
    return (!subjIsReverseComplemented(m, options, TFormat()))
            ? (m.subjEnd   * 3 - getQryFrameShift(m, options, TFormat()))
            : (m.subjStart * 3 + getQryFrameShift(m, options, TFormat()));
}

template <BlastFormatOptions::M mf,
          BlastFormatOptions::Generation g>
constexpr auto
getTrueSubjEndPos(Match const & m,
             ScarletOptions const & options,
             BlastFormat<mf,BlastFormatOptions::TBlastX,g> const & /*tag*/)
-> decltype(m.subjEnd)
{
    typedef BlastFormat<mf,BlastFormatOptions::TBlastX,g> TFormat;
    return (!subjIsReverseComplemented(m, options, TFormat()))
            ? (m.subjEnd   * 3 - getQryFrameShift(m, options, TFormat()))
            : (m.subjStart * 3 + getQryFrameShift(m, options, TFormat()));
}



#endif // header guard