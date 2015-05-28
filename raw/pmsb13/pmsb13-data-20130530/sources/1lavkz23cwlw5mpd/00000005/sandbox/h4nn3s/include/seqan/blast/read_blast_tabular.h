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
// This file contains routines to read BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_READ_BLAST_TABULAR_H_
#define SEQAN_EXTRAS_BLAST_READ_BLAST_TABULAR_H_

/* IMPLEMENTATION NOTES

BLAST Tabular example:

The format of a blast tabular output file is less simple than it looks, here's
the general form

HEADER
 RECORD
 RECORD
 RECORD
 RECORD
HEADER
 RECORD
HEADER
HEADER
...

=> Header for each sequence, 0-n Records for each sequence
=> Each record is one-line, each Header is multiline

A Header usually consists of:

# Program Tag [e.g BLASTX 2.2.27+ ]
# Query Id [ID of the query *sequence*]
# Database Id [note that this is not the name of the sequence in the db, but of
  the database itself, e.g. "nr" -> usually the same for each file]
# Fields: [Headers of columns]
# n "hits found"

The first three lines are always written.
The Fields line is always writen by NCBI Blast, but only when hits > 0 by NCBI Blast+.
The "number of hits"-line is always printed by NCBI Blast+, and never by NCBI Blast.

Possibly other lines can be written as comments.

Because 0 records are allowed, multiple Headers can succeed each other, the criterium for seperation employed by this implementation is that an NCBI Blast
records always ends after the "Fields" line and NCBI Blast+ records end after
the "number of hits"-line.
A file is considered NCBI Blast+ format, when it has comment line that starts
with "# BLAST" and ends with "+". In all other cases it is considered
traditional Blast format.
*/


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

// ============================================================================
// Functions
// ============================================================================v //

//TODO(h4nn3s): doc
template <typename TFile,
          typename TPass,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline bool
onRecord(RecordReader<TFile, TPass > & reader,
         BlastFormat<BlastFormatOptions::Tabular,
                     p,
                     g> const & /*tag*/)
{
    if (value(reader) == '#')
        return false;
    return true;
}

template <typename TFile,
          typename TPass,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline bool
onRecord(RecordReader<TFile, TPass > & reader,
         BlastFormat<BlastFormatOptions::TabularWithHeader,
                     p,
                     g> const & /*tag*/)
{
    return onRecord(reader,
                    BlastFormat<BlastFormatOptions::Tabular,p,g>());
}


// verify whether the fields are default fields
template <typename TString,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_verifyFields(StringSet<TString> const & fields,
              unsigned const hits, // irrelevant for traditional header
              BlastFormat<BlastFormatOptions::TabularWithHeader,
                          p,
                          g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        BlastFormatOptions::Blast>  TFormat;

    if (g == BlastFormatOptions::BlastPlus)
        if ((hits == 0) && length(fields) )
            return 0;


    CharString fieldStr;
    joinStringSet(fieldStr, fields, ", ");

    // it is only relevant that the firtst 12 fields by as expected
    return prefix(fieldStr, length(_defaultFields(TFormat())))
              == _defaultFields(TFormat());
}

// ----------------------------------------------------------------------------
// Function readHeader()                               [Single pass]
// ----------------------------------------------------------------------------

template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFile,
          typename TPass,
          typename TString,
          typename TString2,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_readHeaderImplBlastTab(TqId & qId,
                        TDBName & dbName,
                        TVersionString & versionString,
                        StringSet<TString> & fields,
                        unsigned long & hits,
                        StringSet<TString2> & otherLines, // any other lines
                        RecordReader<TFile, TPass > & reader,
                        const bool strict,
                        BlastFormat<BlastFormatOptions::TabularWithHeader,
                                    p,
                                    g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> TFormat;

    // this is a record instead of a header
    if (onRecord(reader, TFormat()))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    int ret = 0;
    bool keyCorrect = false;
    int queryLinePresent = 0;
    int dbLinePresent = 0;
    int fieldsLinePresent = 0;
    int hitsLinePresent = 0;
    int lastLine = 0;

    while ((!atEnd(reader)) && (!onRecord(reader, TFormat())))// in Header
    {
        // skip '#'
        ret = goNext(reader);
        if (ret)
            return ret;
        //skip blanks following the '#'
        ret = skipBlanks(reader);
        if (ret)
            return ret;

        CharString key = "";
        ret = readUntilWhitespace(key, reader);
        if (ret)
            return ret;

        if (key == _programTagToString(TFormat()))
        {
            // get whole line
            CharString buf;
            ret = readLine(buf, reader);
            if (ret)
                return ret;
            CharString &fullLine = key;
            append(fullLine, buf);

            versionString = fullLine;
            keyCorrect = true;
        }
        else if (key == "Query:")
        {
            ret = skipBlanks(reader);
            if (ret)
                return ret;
            ret = readLine(qId, reader);
            if (ret)
                return ret;
            ++queryLinePresent;
        }
        else if (key == "Database:")
        {
            ret = skipBlanks(reader);
            if (ret)
                return ret;
            ret = readLine(dbName, reader);
            if (ret)
                return ret;
            ++dbLinePresent;
        }
        else if (key == "Fields:")
        {
            ret = skipBlanks(reader);
            if (ret)
                return ret;

            CharString buf;
            ret = readLine(buf, reader);
            if (ret)
                return ret;
            strSplit(fields, buf, ", ");

            ++fieldsLinePresent;
            if (g == BlastFormatOptions::Blast)
                break; // header is finished
        }
        else
        {
            // get whole line
            CharString buf;
            ret = readLine(buf, reader);
            if (ret)
                return ret;
            CharString &fullLine = key;
            append(fullLine, buf);

            if (g == BlastFormatOptions::BlastPlus)
            {
                // last line of BlastPlus Format
                if (hasPrefix(fullLine, "BLAST processed"))
                    ++lastLine;
                // is hits counter?
                else if (hasSuffix(fullLine, "hits found"))
                {
                    CharString hitsString = "";
                    for (unsigned i = 0;
                        (i < length(fullLine) && isdigit(fullLine[i]));
                        ++i)
                        append(hitsString, fullLine[i], Generous());

                    ret = !lexicalCast2(hits, hitsString);
                    if (ret && strict)
                        return RecordReader<TFile, TPass >::INVALID_FORMAT;

                    ++hitsLinePresent;
                    break; // header is finished
                }
            }
            else
                appendValue(otherLines, fullLine, Generous());
        }
    }

    if (!strict)
        return 0;

    if (g == BlastFormatOptions::Blast)
        if (  keyCorrect                &&
             (queryLinePresent   == 1)  &&
             (dbLinePresent      == 1)  &&
             (fieldsLinePresent  == 1)  &&
             (length(otherLines) == 0)   )
            return 0;

    if (g == BlastFormatOptions::BlastPlus)
        if (( keyCorrect                                   &&
             (queryLinePresent    == 1)                    &&
             (dbLinePresent       == 1)                    &&
             (hitsLinePresent     == 1)                    &&
             ( (fieldsLinePresent == 1) || (hits==0) )     &&
             (length(otherLines)  == 0)                     ) ||
             ( (lastLine          == 1) && atEnd(reader))   )
            return 0;

    return RecordReader<TFile, TPass >::INVALID_FORMAT;
}

/**
.Function.BLAST#readHeader
..cat:Input / Output
..signature:readHeader(qId, dbName, versionString, [hits,] [fields, [otherLines,]] recordReader, strict, BlastFormat)
..param.qId: String to hold the query ID from the header
...type:Class.String
...type:nolink:or similar
..param.dbname: String to hold the database name from the header
...type:Class.String
...type:nolink:or similar
..param.versionString: String to hold the Blast program Tag and Version
...type:Class.String
...type:nolink:or similar
..param.hits: Numerical to hold the number of hits that will follow the header (only available in BlastPlus spec of BlastFormat)
...type:nolink:$unsigned long$
..param.fields: StringSet to hold column identifiers, useful if non-defaults are expected
...type:Class.StringSet
..param.otherLines: StringSet to hold any comment or header lines that are not identified otherwise
...type:Class.StringSet
..param.recodReader: the RecordReader we are reading from
...type:Class.RecordReader
..param.strict: switch to signify whether the function should return error on a non-conforming header or just "get whatever possible". If not using strict, it is recommended to pass fields and otherLines and verify these manually.
..type:nolink:$bool$
..param.BlastFormat: The format tag, note that BlastFormat must be further specified
...type:Class.BlastFormat
..summary: read a Header from a Blast output file
..remark: call this function whenever a line in the recordReader begins with '#'
..include:seqan/blast.h
*/

// default traditional Blast or BlastPlus
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFile,
          typename TPass,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
readHeader(TqId & qId,
           TDBName & dbName,
           TVersionString & versionString,
           RecordReader<TFile, TPass > & reader,
           const bool strict,
           BlastFormat<BlastFormatOptions::TabularWithHeader,
                       p,
                       g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> TFormat;

    StringSet<CharString> otherLines;
    StringSet<CharString> fields;
    unsigned long hits = 0;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    if (strict)
    {
        ret =  _verifyFields(fields, hits, TFormat());
        if (ret)
            RecordReader<TFile, TPass >::INVALID_FORMAT;
    }

    return 0;
}

// BlastPlus with hit count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFile,
          typename TPass,
          BlastFormatOptions::Program p>
inline int
readHeader(TqId & qId,
           TDBName & dbName,
           TVersionString & versionString,
           unsigned long & hits,
           RecordReader<TFile, TPass > & reader,
           const bool strict,
           BlastFormat<BlastFormatOptions::TabularWithHeader,
                       p,
                       BlastFormatOptions::BlastPlus> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        BlastFormatOptions::BlastPlus> TFormat;

    StringSet<CharString> otherLines;
    StringSet<CharString> fields;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    if (strict)
    {
        ret =  _verifyFields(fields, hits, TFormat());
        if (ret)
            RecordReader<TFile, TPass >::INVALID_FORMAT;
    }

    return 0;
}


// with fields
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldString,
          typename TFile,
          typename TPass,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
readHeader(TqId & qId,
           TDBName & dbName,
           TVersionString & versionString,
           StringSet<TFieldString> & fields,
           RecordReader<TFile, TPass > & reader,
           const bool strict,
           BlastFormat<BlastFormatOptions::TabularWithHeader,
                       p,
                       g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> TFormat;

    StringSet<CharString> otherLines;
    unsigned long hits =0;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
    return 0;
}

// with fields and hits count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldString,
          typename TFile,
          typename TPass,
          BlastFormatOptions::Program p>
inline int
readHeader(TqId & qId,
           TDBName & dbName,
           TVersionString & versionString,
           unsigned long & hits,
           StringSet<TFieldString> & fields,
           RecordReader<TFile, TPass > & reader,
           const bool strict,
           BlastFormat<BlastFormatOptions::TabularWithHeader,
                       p,
                       BlastFormatOptions::BlastPlus> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        BlastFormatOptions::BlastPlus> TFormat;

    StringSet<CharString> otherLines;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
    return 0;
}


// with fields and otherLines
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldString,
          typename TOtherString,
          typename TFile,
          typename TPass,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
readHeader(TqId & qId,
           TDBName & dbName,
           TVersionString & versionString,
           StringSet<TFieldString> & fields,
           StringSet<TOtherString> & otherLines,
           RecordReader<TFile, TPass > & reader,
           const bool strict,
           BlastFormat<BlastFormatOptions::TabularWithHeader,
                       p,
                       g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> TFormat;

    unsigned long hits =0;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
    return 0;
}

// with fields and otherLines and hits count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldString,
          typename TOtherString,
          typename TFile,
          typename TPass,
          BlastFormatOptions::Program p>
inline int
readHeader(TqId & qId,
           TDBName & dbName,
           TVersionString & versionString,
           unsigned long & hits,
           StringSet<TFieldString> & fields,
           StringSet<TOtherString> & otherLines,
           RecordReader<TFile, TPass > & reader,
           const bool strict,
           BlastFormat<BlastFormatOptions::TabularWithHeader,
                       p,
                       BlastFormatOptions::BlastPlus> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        BlastFormatOptions::BlastPlus> TFormat;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
    return 0;
}

// ----------------------------------------------------------------------------
// Function skipHeader()                               [Single pass]
// ----------------------------------------------------------------------------

/**
.Function.BLAST#skipHeader
..cat:Input / Output
..signature:skipHeader(recordReader, [strict,] BlastFormat)
..param.recodReader: the RecordReader we are reading from
...type:Class.RecordReader
..param.strict: switch to signify whether the function should return error on a non-conforming header or just "get whatever possible". False if not specified.
..type:nolink:$bool$
..param.BlastFormat: The format tag, note that BlastFormat must be further specified
...type:Class.BlastFormat
..summary: skip a header from a Blast tabular output file, optionally verifying it for format compliance.
..remark: call this function whenever you want to skip exactly one header. If you want to go directly to the beginning of the next record (possibly skipping multiple headers that have no succeeding records) use @Function.BLAST#skipUntilRecord@ instead.
..include:seqan/blast.h
*/
template <typename TFile,
          typename TPass,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
skipHeader(RecordReader<TFile, TPass > & reader,
           const bool strict,
           BlastFormat<m,p,g> const & /*tag*/)
{
    CharString qId;
    CharString dbName;
    CharString versionString;
    StringSet<CharString> fields;
    unsigned long hits = 0;
    StringSet<CharString> otherLines;

    return  _readHeaderImplBlastTab(qId,
                                    dbName,
                                    versionString,
                                    fields,
                                    hits,
                                    otherLines,
                                    reader,
                                    strict,
                                    BlastFormat<m,p,g>());
}

template <typename TFile,
          typename TPass,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
skipHeader(RecordReader<TFile, TPass > & reader,
           BlastFormat<m,p,g> const & /*tag*/)
{
    return skipHeader(reader, false, BlastFormat<m,p,g>());
}


// ----------------------------------------------------------------------------
// Function skipUntilRecord()                           [Single pass]
// ----------------------------------------------------------------------------

/**
.Function.BLAST#skipUntilRecord
..cat:Input / Output
..signature:skipUntilRecord(recordReader, BlastFormat)
..param.recodReader: the RecordReader we are reading from
...type:Class.RecordReader
..param.BlastFormat: The format tag, note that BlastFormat must be further specified
...type:Class.BlastFormat
..summary: skip arbitrary number of headers and/or comment lines until a record is reached
..remark: call this function whenever you are on a comment character ('#') in the file and want to jump to the beginning of the next record. If you want skip only a single header (to count skipped headers or verify its conformance
to standards), use @Function.BLAST#skipHeader@ instead.
..include:seqan/blast.h
*/
template <typename TFile,
          typename TPass,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
skipUntilRecord(RecordReader<TFile, TPass > & reader,
                BlastFormat<m,p,g> const & /*tag*/)
{
    int ret = 0;
    while ((!atEnd(reader)) && value(reader) == '#') // skip comments
    {
        ret = skipLine(reader);
        if (ret)
            return ret;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                               [Single pass]
// ----------------------------------------------------------------------------

template <typename TString,
          typename TFile,
          typename TPass>
inline int
_readRecordImplBlastTabArbitrary(StringSet<TString> & fields,
                                 RecordReader<TFile, TPass > & reader)
{
    CharString buf;
    int ret = readLine(buf, reader);
    if ( (ret) && (!atEnd(reader)) ) // file could end w/o empty newline
        return ret;

    strSplit(fields, buf, '\t', false);
    return 0;
}

template <typename TqId,
          typename TsId,
          typename TFile,
          typename TPass>
inline int
_readRecordImplBlastTabDefault(TqId & qId,
                               TsId & sId,
                               double  & percentIdent,
                               unsigned & ali_length,
                               unsigned & num_mismatches,
                               unsigned & gap_openings,
                               unsigned long & qStart,
                               unsigned long & qEnd,
                               unsigned long & sStart,
                               unsigned long & sEnd,
                               double & eval,
                               double & bitScore,
                               RecordReader<TFile, TPass > & reader)
{
    int ret = 0;

    for (int i = 0; i < 12; ++i)
    {
        CharString buf;
        if (i < 11)
        {
            ret = readUntilChar(buf, reader, '\t');
            if (ret)
                return ret;
            ret = goNext(reader); //skip '\t', go to begin of next field
            if (ret)
                return ret;
        } else
        {
            ret = readUntilTabOrLineBreak(buf, reader);
            if (ret)
                return ret;
            ret = skipLine(reader); // skip extra fields or '\n'
            if ( (ret) && (!atEnd(reader)) ) // file may end w/o empty newline
                return ret;
        }

        switch(i)
        {
            case  0: ret = lexicalCast2(qId, buf); break;
            case  1: ret = lexicalCast2(sId, buf); break;
            case  2: ret = lexicalCast2(percentIdent, buf); break;
            case  3: ret = lexicalCast2(ali_length, buf); break;
            case  4: ret = lexicalCast2(num_mismatches, buf); break;
            case  5: ret = lexicalCast2(gap_openings, buf); break;
            case  6: ret = lexicalCast2(qStart, buf); break;
            case  7: ret = lexicalCast2(qEnd, buf); break;
            case  8: ret = lexicalCast2(sStart, buf); break;
            case  9: ret = lexicalCast2(sEnd, buf); break;
            case 10: ret = lexicalCast2(eval, buf); break;
            case 11: ret = lexicalCast2(bitScore, buf); break;
        }
        ret = !ret; // lexicalCast2 returns true when successful
        if (ret)
            return RecordReader<TFile, TPass >::INVALID_FORMAT;
    }
    return 0;
}


/**
.Function.BLAST#readRecord
..cat:Input / Output
..signature:readRecord(qId, sId, percentIdent, aliLength, numMismatches, gapOpenings, qStart, qEnd, sStart, sEnd, eValue, bitScore, recordReader, BlastFormat)
..signature:readRecord(fields, recordReader, BlastFormat)
..param.qId: ID of the query
...type:Class.String
...type:nolink:or similar
..param.sId: ID of the subject (sequence in database)
...type:Class.String
...type:nolink:or similar
..param.percentIdent: percentage of identies in alignment
...type:nolink:double
..param.aliLength: length of alignment
...type:nolink:unsigned
..param.numMismatches: number of mismatches in alignment
...type:nolink:unsigned
..param.gapOpenings: number consecutive gaps per alignment
...type:nolink:unsigned
..param.qStart: alignment begin position on query
...type:nolink:unsigned long
..param.qEnd: alignment end position on query
...type:nolink:unsigned long
..param.sStart: alignment begin position on subject
...type:nolink:unsigned long
..param.sEnd: alignment end position on subject
...type:nolink:unsigned long
..param.eValue: alignment e-Value
...type:nolink:double
..param.bitScore: alignment bit-Score
...type:nolink:double
..param.recodReader: the RecordReader we are reading from
...type:Class.RecordReader
..param.BlastFormat: The format tag, note that BlastFormat must be further specified
...type:Class.BlastFormat
..summary: skip arbitrary number of headers and/or comment lines until a record is reached
..remark: call this function whenever you are on a comment character ('#') in the file and want to jump to the beginning of the next record. If you want skip only a single header (to count skipped headers or verify its conformance
to standards), use @Function.BLAST#skipHeader@ instead.
..include:seqan/blast.h
*/
template <typename TqId,
          typename TsId,
          typename TFile,
          typename TPass,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
readRecord(TqId & qId,
           TsId & sId,
           double  & percentIdent,
           unsigned & ali_length,
           unsigned & num_mismatches,
           unsigned & gap_openings,
           unsigned long & qStart,
           unsigned long & qEnd,
           unsigned long & sStart,
           unsigned long & sEnd,
           double & eval,
           double & bitScore,
           RecordReader<TFile, TPass > & reader,
           BlastFormat<m,p,g> const & /*tag*/)
{
    // header should have been read or skipped
    if (!onRecord(reader, BlastFormat<m,p,g>()))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    return _readRecordImplBlastTabDefault(qId,
                                          sId,
                                          percentIdent,
                                          ali_length,
                                          num_mismatches,
                                          gap_openings,
                                          qStart,
                                          qEnd,
                                          sStart,
                                          sEnd,
                                          eval,
                                          bitScore,
                                          reader);
}

template <typename TFile,
          typename TPass,
          typename TString,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
readRecord(StringSet<TString> & fields,
           RecordReader<TFile, TPass > & reader,
           BlastFormat<m,p,g> const & /*tag*/)
{
    // header should have been read or skipped
    if (!onRecord(reader, BlastFormat<m,p,g>()))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    return _readRecordImplBlastTabArbitrary(fields, reader);
}

} // namespace seqan
#endif // header guard
