// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Module for handling NCBI Blast I/O and E-Value computation
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_BLAST_BASE_H_
#define SEQAN_EXTRAS_BLAST_BLAST_BASE_H_

#if __cplusplus <= 199711L //C++98
#define constexpr inline
#endif

namespace seqan {

template <typename TLeft, typename TRight >
inline bool
hasSuffix(TLeft const & left,
          TRight const & right)
{
    if (length(right) > length(left))
        return false;

    return isEqual(suffix(left, length(left)-length(right)), right);

}

template <typename TLeft, typename TRight >
inline bool
isSuffix(TLeft const & left,
          TRight const & right)
{
    return hasSuffix(right, left);
}

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BlastFormatOptions
{
/**
.Enum.BlastFormatOptions\colon\colonM
..cat:Input / Output
..summary:Enum with BLAST file format specs
..value.Pairwise: default blast output (blastall -m 0 & blast* -outfmt 0)
..value.MasterSlaveIdent: master-slave showing identities (blastall -m 1 & blast* -outfmt 1)
..value.MasterSlaveNoIdent: master-slave without identities (blastall -m 2 & blast* -outfmt 2)
..value.FlatMasterSlaveIdent: flat master-slave showing identities (blastall -m 3 & blast* -outfmt 3)
..value.FlatMasterSlaveNoIdent: flat master-slave without identities (blastall -m 4 & blast* -outfmt 4)
..value.MasterSlaveBluntEnds: master-slave without identities, with blunt ends (blastall -m 5, not available with Blast+)
..value.FlatMasterSlaveBluntEnds: flat master-slave without identities, with blunt ends (blastall -m 6, not available with Blast+)
..value.XML: XML (blastall -m 7 & blast* -outfmt 8)
..value.Tabular: tab-seperated (blastall -m 8 & blast* -outfmt 6)
..value.TabularWithHeader: tab-seperated with Header / comments (blastall -m 9 & blast* -outfmt 7)
..value.TextASN1: Abstract Syntax Notation One (blast* -outfmt 8, not available in traditional BLAST)
..value.BinASN1: Abstract Syntax Notation One (blast* -outfmt 9, not available in traditional BLAST)
..value.CSV: comma-seperated values (blast* -outfmt 10, not available in traditional BLAST)
..value.BlastArASN1: Blast Archive Format, Abstract Syntax Notation One (blast* -outfmt 11, not available in traditional BLAST)
..remark: SeqAn currently implements Pairwise(TODO!), Tabular and TabularWithHeader
..include:seqan/blast.h
*/
    enum M
    {
        Pairwise = 0,
        MasterSlaveIdent = 1,
        MasterSlaveNoIdent = 2,
        FlatMasterSlaveIdent = 3,
        FlatMasterSlaveNoIdent = 4,
        MasterSlaveBluntEnds = 5,       // only available in Generation==Blast
        FlatMasterSlaveBluntEnds = 6,   // only available in Generation==Blast
        XML = 7,
        Tabular = 8,
        TabularWithHeader = 9,
        TextASN1 = 10,                 // only available in Generation==Blast+
        BinASN1 = 11,                  // only available in Generation==Blast+
        CSV = 12,                      // only available in Generation==Blast+
        BlastArASN1 = 13,              // only available in Generation==Blast+
        INVALID_M=1023
    };
/**
.Enum.BlastFormatOptions\colon\colonProgram
..cat:Input / Output
..summary:Enum with BLAST program spec
..value.BlastN: Nucleotide Query VS Nucleotide Subject
..value.BlastP: Protein Query VS Protein Subject
..value.BlastX: translated Nucleotide Query VS Protein Subject
..value.TBlastN: Protein Query VS translated Nucleotide Subject
..value.TBlastX: translated Nucleotide Query VS translated Nucleotide Subject
..include:seqan/blast.h
*/
    enum Program
    {
        BlastN,         //              nucl VS             nucl
        BlastP,         //              prot VS             prot
        BlastX,         // translated   nucl VS             prot
        TBlastN,        //              prot VS translated  nucl
        TBlastX,        // translated   nucl VS translated  nucl
        INVALID_Program=1023
    };
/**
.Enum.BlastFormatOptions\colon\colonGeneration
..cat:Input / Output
..summary:Enum with BLAST program Generation
..value.Blast: traditional NCBI Blast, written in C
..value.BlastPlus: NCBI Blast+, written in C++
..include:seqan/blast.h
*/
    enum Generation
    {
        Blast,
        BlastPlus,
        INVALID_Generation=1023
    };
};

/**
.Class.BlastFormat
..cat:Input / Output
..summary:Blast Format specifier
..signature:BlastFormat<m, p, g>
..param.m: File Type Format
...type:Enum.BlastFormatOptions::M
..param.p: Program Type Format
...type:Enum.BlastFormatOptions::Program
..param.g: Program Generation
...type:Enum.BlastFormatOptions::Generation
..include:seqan/blast.h
*/
/*TODO(C++11): change struct BlastFormat to struct BlastFormat_ and wrap
a type-dependent typedef Tag around it */
template <BlastFormatOptions::M            _m,
          BlastFormatOptions::Program      _p,
          BlastFormatOptions::Generation   _g>
struct BlastFormat
{
    // have static members for run-time acces to "type"
    static const BlastFormatOptions::M          m = _m;
    static const BlastFormatOptions::Program    p = _p;
    static const BlastFormatOptions::Generation g = _g;

};

// ============================================================================
// Metafunctions
// ============================================================================
/*
//TODO doc
template<BlastFormatOptions::Program p = BlastFormatOptions::INVALID_Program>
struct BlastInputType
{
    enum Alph
    {
        Nucl,
        Prot,
        TransNucl
    };

};

template<>
struct BlastInputType<BlastFormatOptions::BlastN>
{
    enum Alpha
    {
        QUERY = BlastInputType::Nucl,
        SUBJ = BlastInputType::Nucl
    };
};

template<>
struct BlastInputType<BlastFormatOptions::BlastP>
{
    enum Alpha
    {
        QUERY = BlastInputType::Prot,
        SUBJ = BlastInputType::Prot
    };
};

template<>
struct BlastInputType<BlastFormatOptions::BlastX>
{
    enum Alpha
    {
        QUERY = BlastInputType::TransNucl,
        SUBJ = BlastInputType::Prot
    };
};

template<>
struct BlastInputType<BlastFormatOptions::TBlastN>
{
    enum Alpha
    {
        QUERY = BlastInputType::Prot,
        SUBJ = BlastInputType::TransNucl
    };
};

template<>
struct BlastInputType<BlastFormatOptions::TBlastX>
{
    enum Alpha
    {
        QUERY = BlastInputType::TransNucl,
        SUBJ = BlastInputType::TransNucl
    };
};
*/
// ----------------------------------------------------------------------------
// getBlastProgramType()
// ----------------------------------------------------------------------------

//TODO(h4nn3s): test, document
template< typename TQueryAlph, typename TSubjAlph>
constexpr
BlastFormatOptions::Program
getBlastProgramType(const TQueryAlph &, const TSubjAlph &)
{
    return BlastFormatOptions::INVALID_Program;
}

template<typename TQueryAlph, typename TSubjAlph,
         typename TSpec, typename TSpec2>
constexpr
BlastFormatOptions::Program
getBlastProgramType(const String<TQueryAlph, TSpec> &,
                    const String<TSubjAlph, TSpec2> &)
{
    // needs constexpr constructors of Alphabet types
    return getBlastProgramType(TQueryAlph(), TSubjAlph());
}

// --- DNA vs DNA ---
// NOTE that Dna VS Dna could also be TBlastX, but BlastN is more common
constexpr
BlastFormatOptions::Program
getBlastProgramType(const Dna &, const Dna &)
{
    return BlastFormatOptions::BlastN;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(const Dna &, const Dna5 &)
{
    return BlastFormatOptions::BlastN;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(const Dna5 &, const Dna &)
{
    return BlastFormatOptions::BlastN;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(const Dna5 &, const Dna5 &)
{
    return BlastFormatOptions::BlastN;
}

// --- Protein vs Protein ---
constexpr
BlastFormatOptions::Program
getBlastProgramType(const AminoAcid &, const AminoAcid &)
{
    return BlastFormatOptions::BlastP;
}

// --- Dna vs Protein ---
constexpr
BlastFormatOptions::Program
getBlastProgramType(const Dna &, const AminoAcid &)
{
    return BlastFormatOptions::BlastX;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(const Dna5 &, const AminoAcid &)
{
    return BlastFormatOptions::BlastX;
}

// --- Protein vs Dna ---
constexpr
BlastFormatOptions::Program
getBlastProgramType(const AminoAcid &, const Dna &)
{
    return BlastFormatOptions::TBlastX;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(const AminoAcid &, const Dna5 &)
{
    return BlastFormatOptions::TBlastX;
}




// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// _programTagToString()
// ----------------------------------------------------------------------------


//TODO(h4nn3s): document?

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             p,
                                             g> const & /*tag*/)
{
    return "UNKOWN BLAST PROGRAM";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::BlastN,
                                             g> const & /*tag*/)
{
    return "BlastN";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::BlastP,
                                             g> const & /*tag*/)
{
    return "BLASTP";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::BlastX,
                                             g> const & /*tag*/)
{
    return "BLASTX";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::TBlastN,
                                             g> const & /*tag*/)
{
    return "TBLASTN";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::TBlastX,
                                             g> const & /*tag*/)
{
    return "TBLASTX";
}

// ----------------------------------------------------------------------------
// _seqanProgramTag()
// ----------------------------------------------------------------------------


template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline
const char * _seqanProgramTag(BlastFormat<m,p,g> const & tag)
{
    // get quotes around numerical MACROS
    #define QUOTE_X(t)#t
    #define QUOTE(t)QUOTE_X(t)

    // make nice, long string at compile time
    #define STR " I/O Module of SeqAn-" QUOTE(SEQAN_VERSION_MAJOR) "." QUOTE(SEQAN_VERSION_MINOR) "." QUOTE(SEQAN_VERSION_PATCH) " (http://www.seqan.de)"

    CharString s = _programTagToString(tag);
    append(s, STR);
    return toCString(s);

    //TODO(C++): SeqAn could implement some proper compile-time string-handling
    // which would enable this function to be constexpr

    #undef STR
    #undef QUOTE
    #undef QUOTE_X
}

// ----------------------------------------------------------------------------
// _defaultFields()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::Generation g>
constexpr
const char * _defaultFields()
{
    return "ERROR Fields not specializied for this type";
}

template <>
constexpr
const char * _defaultFields<BlastFormatOptions::Blast>()
{
    return "Fields: Query id, Subject id, % identity, alignment length," \
           " mismatches, gap openings, q. start, q. end, s. start, s." \
           " end, e-value, bit score";
}

template <>
constexpr
const char * _defaultFields<BlastFormatOptions::BlastPlus>()
{
    return "Fields: query id, subject id, % identity, alignment " \
           "length, mismatches, gap opens, q. start, q. end, s. " \
           "start, s. end, evalue, bit score";
}

template <BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
constexpr
const char * _defaultFields(BlastFormat<BlastFormatOptions::TabularWithHeader,
                                        p,
                                        g> const &)
{
    return _defaultFields<g>();
}


} // namespace seqan

#if __cplusplus <= 199711L //C++98
#undef constexpr
#endif

#endif // header guard
