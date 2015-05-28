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
// This file contains routines to generate BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_WRITE_BLAST_TABULAR_H_
#define SEQAN_EXTRAS_BLAST_WRITE_BLAST_TABULAR_H_

#include <sstream>
#include <seqan/version.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

#ifdef SEQAN_CPP11
template <BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
constexpr
const char *
_seperatorString(BlastFormat<BlastFormatOptions::TabularWithHeader,
                              p,
                              g> const & /*tag*/)
{
    return ", ";
}

template <BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
constexpr
const char *
_seperatorString(BlastFormat<BlastFormatOptions::Tabular,
                              p,
                              g> const & /*tag*/)
{
    return "\t";
}


#endif //SEQAN_CPP11

// ============================================================================
// Functions
// ============================================================================

#ifdef SEQAN_CPP11
// ----------------------------------------------------------------------------
// Helper functions for printing n columns (requires C++11 variadic templates)
// ----------------------------------------------------------------------------


template <typename TStream, typename TTag>
inline int
_writeFields(TStream & stream, TTag const & /*tag*/)
{
    return streamPut(stream, '\n');
}

template <typename TStream, typename TField, typename... TFields, typename TTag>
inline int
_writeFields(TStream & stream,
             TTag const & /*tag*/,
             TField const & field1, const TFields&... fields)
{
    int ret = streamPut(stream,  _seperatorString(TTag()));
    if (ret)
        return ret;

    ret = streamPut(stream, field1);
    if (ret)
        return ret;

    return _writeFields(stream, TTag(), fields... );
}
#endif //SEQAN_CPP11

// ----------------------------------------------------------------------------
// Function writeHeader()
// ----------------------------------------------------------------------------

/**
.Function.BLAST I/O#writeHeader
..signature:int writeHeader(stream, query_id, db_name, [fields,] BlastFormat)
..signature:int writeHeader(stream, query_id, db_name, BlastFormat[, field1, ... fieldN])
..param.stream:The stream to write to.
...type:Concept.Stream
..param.query_id:ID of the query sequence.
...type:nolink:String-type
..param.db_name:Name of the database / file name.
...type:nolink:String-type
..param.fields:A StringSet with column identifiers to print
...type;StringSet
..param.BlastFormat: The format tag, note that BlastFormat must be further specified
...type:Class.BlastFormat
..param.fieldN:Manually supply headers of columns, defaults to Blast's defaul 12 columns (only supported with C++11)
...type:nolink:String-type
..remarks:For BlastFormat-types that have no header, this is a NOOP.
..include:seqan/blast.h
*/

template <typename TStream, typename TString,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeHeader(TStream &, TString const &, TString const &,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    return 0;
}


template <typename TStream, typename TString,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
_writeHeaderWithoutFields(TStream & stream,
                          TString const & qryId, TString const & dbName,
                          BlastFormat<BlastFormatOptions::TabularWithHeader,
                                      p,
                                      g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,p,g> TFormat;

    std::stringstream ss (std::stringstream::in | std::stringstream::out);

    ss << "# " << _seqanProgramTag(TFormat()) << '\n'
       << "# Query: " << qryId << '\n'
       << "# Database: " << dbName << '\n';

    return streamPut(stream, ss);
}

// default case
template <typename TStream, typename TString,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeHeader(TStream & stream,
            TString const & qryId, TString const & dbName,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,p,g> TFormat;

    int ret = _writeHeaderWithoutFields(stream, qryId, dbName, TFormat());
    if (ret)
        return ret;

    streamPut(stream, "# ");
    if (ret)
        return ret;

    ret = streamPut(stream, _defaultFields(TFormat()));
    if (ret)
        return ret;

    return streamPut(stream, '\n');
}

// 
template <typename TStream, typename TqId, typename TdbName,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeHeader(TStream & stream,
            TqId const & qryId, TdbName const & dbName,
            StringSet<CharString> const & fields,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,p,g> TFormat;

    int ret = _writeHeaderWithoutFields(stream, qryId, dbName, TFormat());
    if (ret)
        return ret;

    std::stringstream ss (std::stringstream::in | std::stringstream::out);

    ss << "# Fields: ";
    for (unsigned i = 0; i < length(fields); ++i)
        ss << fields[i] << ((i == length(fields) -1)
                                ? "\n"
                                : _seperatorString(TFormat()));

    return streamPut(stream, ss);
}

#ifdef SEQAN_CPP11
// Functions for arbitrary number and typed fields

template <typename TStream, typename TString,
          typename TField, typename... TFields,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeHeader(TStream & stream,
            TString const & qryId, TString const & dbName,
            BlastFormat<m,p,g> const & /*tag*/,
            TField const & field1,
            const TFields&... fields)
{
    return 0;
}

template <typename TStream, typename TString,
          typename TField, typename... TFields
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeHeader(TStream & stream,
            TString const & qryId, TString const & dbName,
            BlastFormat<m,p,g> const & /*tag*/,
            TField const & field1,
            const TFields&... fields)
{
    typedef BlastFormat<m,p,g> Format;

    int ret = _writeHeaderWithoutFields(stream, qryId, dbName, Format());
    if (ret)
        return ret;

    ret = streamPut(stream, "# Fields: ");
    if (ret)
        return ret;

    ret = streamPut(stream, field1);

    return _writeFields(stream, Format(), fields...);
}

#endif //SEQAN_CPP11

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/**
.Function.BLAST I/O#writeRecord
..signature:int writeRecord(stream, query_id, subject_id, percentIdent, ali_length, num_mismatches, gap_openings, qStart, qEnd, sStart, sEnd, eval, bitScore, BlastFormat)
..signature:int writeRecord(stream, fields, BlastFormat)
..signature:int writeRecord(stream, query_id, subject_id, fields, BlastFormat)
..signature:int writeRecord(stream, BlastFormat, field1,[ ... fieldN])
..param.stream:The stream to write to.
...type:Concept.Stream
..param.query_id:ID of the query sequence.
...type:nolink:String-type
..param.subject_id:ID of the database sequence.
...type:nolink:String-type
..param.num_identies:number of identical position in the alignment
...type:nolink:double
..param.ali_length:length of alignment
...type:nolink:unsigned
..param.num_mismatches: number of mismatched positions in alignment
...type:nolink:unsigned
..param.gap_openings: number of gap_openings in alignment
...type:nolink:unsigned
..param.qStart: begin position of the alignment on the query
...type:nolink:unsigned
..param.qEnd: end position of the alignment on the query
...type:nolink:unsigned
..param.sStart: begin position of the alignment on the subject
...type:nolink:unsigned
..param.sEnd: end position of the alignment on the subject
...type:nolink:unsigned
..param.fields:A String with fields to print (note that these have to be of
the same type, e.g. String(Set) of Strings, or String of Double)
..type:String
..param.BlastFormat: The format tag, note that BlastFormat must be further specified
...type:Class.BlastFormat
..param.fieldN:parameters of differently typed values to print (only supported with C++11)
...type:nolink:anything printable by @streamPut@
..include:seqan/blast.h
*/

template <typename TStream, typename TqId, typename TsId,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeRecord(TStream & stream,
            TqId const & qId, TsId const & sId, unsigned const & num_identies,
            unsigned const & ali_length, unsigned const & num_mismatches,
            unsigned const & gap_openings,
            unsigned long const & qStart, unsigned long const & qEnd,
            unsigned long const & sStart, unsigned long const & sEnd,
            double const & eval, double const & bitScore,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << qId << '\t' << sId << '\t'
       << double(num_identies) * 100 / ali_length << '\t'
       << ali_length << '\t' << num_mismatches << '\t' << gap_openings << '\t'
       << qStart << '\t' << qEnd << '\t' << sStart << '\t' << sEnd << '\t'
       << eval << '\t' << bitScore << '\n';

    return streamPut(stream, ss);
}

template <typename TStream, typename TqId, typename TsId, typename TField,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeRecord(TStream & stream,
            TqId const & qId, TsId const & sId,
            String<TField> const & fields,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    std::stringstream ss (std::stringstream::in | std::stringstream::out);

    ss << qId << '\t' << sId;
    for (int i = 0; i < length(fields); ++i)
        ss << '\t' << fields[i];
    ss << '\n';

    return streamPut(stream, ss);
}

template <typename TStream, typename TField,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeRecord(TStream & stream,
            String<TField> const & fields,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeRecord(stream, "", "", fields, TFormat());
}


// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TqId, typename TsId,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeRecord(TStream & stream,
            TqId const & qId, TsId const & sId, double const & percentIdent,
            unsigned const & ali_length, unsigned const & num_mismatches,
            unsigned const & gap_openings,
            unsigned long const & qStart, unsigned long const & qEnd,
            unsigned long const & sStart, unsigned long const & sEnd,
            double const & eval, double const & bitScore,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeRecord(stream, qId, sId, percentIdent, ali_length,
                       num_mismatches, gap_openings, qStart, qEnd, sStart,
                       sEnd, eval, bitScore, TFormat());
}

// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TqId, typename TsId, typename TField,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeRecord(TStream & stream,
            TqId const & qId, TsId const & sId,
            String<TField> const & fields,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeRecord(stream, qId, sId, fields, TFormat());
}

// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TField,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeRecord(TStream & stream,
            String<TField> const & fields,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeRecord(stream, "", "", fields, TFormat());
}


#ifdef SEQAN_CPP11
// Functions for arbitrary number and typed fields

template <typename TStream, typename TField, typename... TFields,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeRecord(TStream & stream,
            BlastFormat<BlastFormatOptions::Tabular,
                        p,
                        g> const & /*tag*/)
            TField const & field1,
            const TFields&... fields)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    int ret = streamPut(stream, field1);
    if (ret)
        return ret;

    return _writeFields(stream, TFormat(), fields...);
}

// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TField, typename... TFields,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeRecord(TStream & stream,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/,
            TField const & field1,
            const TFields&... fields)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeRecord(stream, TFormat(), field1, fields... );
}
#endif //SEQAN_CPP11




} // namespace seqan
#endif // header guard
