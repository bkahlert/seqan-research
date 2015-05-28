// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the Writer specialization for BLAST Format
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_WRITER_BLAST_H_
#define SEQAN_EXTRAS_MASAI_WRITER_BLAST_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include "tags.h"
#include "store.h"
#include "matches.h"
#include "stream.h"

#include "blast.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


template <typename TScore, typename TProgram,
          typename TSpec=void>
struct WriterBlastSpec
{
    typedef TProgram TProgramType;
    TScore   score;
    typename BlastStatistics<TScore>::TParams    params;

    WriterBlastSpec()
    {}
};


// ----------------------------------------------------------------------------
// Class Writer
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDistance,
          typename TBlastFormat, typename TWriterBlastSpec>
struct Writer<TGenome, TReads, TBlastFormat, TDistance, TWriterBlastSpec>
{
    typedef TBlastFormat                            TFormat;
    typedef Stream<FileStream<char, MMapWriter> >   TStream;
//     typedef ::std::fstream         TStream;

    TGenome                 & genome;
    TReads                  * reads;
    TFragmentStore          & _store;
    TStream                 _stream;
    bool                    disabled;

    TWriterBlastSpec        blastSpec;

    Writer(TGenome & genome, bool disabled = false) :
        genome(genome),
        reads(),
        _store(genome._store),
        disabled(disabled)
    {
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _fillAlignedRead()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function open()                                                     [Writer]
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function close()                                                    [Writer]
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Function setReads()                                                 [Writer]
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function writeAlignments()                                          [Writer]
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Function _resize()                                                  [Writer]
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function _clear()                                                   [Writer]
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Function _writeHeader()                                             [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec>
void _writeHeader(Writer<TGenome, TReads, BlastTab, TDistance, TSpec> & /* writer */)
{
    //TODO(h4nn3s): implement
}

// ----------------------------------------------------------------------------
// Function onMatch()                                          [Writer<Blast>]
// ----------------------------------------------------------------------------


// Single-End
template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Writer<TGenome, TReads, BlastTab, TDistance, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    if (writer.disabled)
        return;

    typedef Writer<TGenome, TReads, BlastTab, TDistance, TSpec>  TWriter;
    typedef Align<TFragmentStore::TReadSeq, ArrayGaps>  TAlign;

    TAlignedReadStoreElement        alignedRead;
    TAlignQualityStoreElement       alignQuality;
    TAlign                          align;

    // Fill aligned read object
    _fillAlignedRead(alignedRead, alignQuality,
                     contigId, beginPos, endPos, readId, errors, reverseComplemented);


    // create alignment
    _alignRead(writer, align, alignedRead, alignQuality, reverseComplemented);

    // construct record
    long sc = 0;
    unsigned int ali_length = 0;
    unsigned int identities = 0;
    double percentIdent = 0;
    unsigned int positives = 0;
    unsigned int mismatches = 0;
    unsigned int gaps = 0;
    unsigned int gap_openings = 0;

    int ret = calcStatsAndScore(sc, ali_length, identities, percentIdent,
                                positives,mismatches, gaps,
                                gap_openings, row(align, 0), row(align, 1),
                                writer.blastSpec.score);
    if (ret)
    {
        ::std::cerr << "Error calculating BLAST stats and score\n";
        return;
    }

    double eval = 0;
    double bitScore = 0;
    const  unsigned long sbjctLength =
                  length(writer._store.contigStore[alignedRead.contigId].seq);
    const unsigned long qryLength =
                  length(writer._store.readSeqStore[alignedRead.readId]);
    ret = calcBitScoreAndEValue(bitScore, eval,
                                sc, qryLength, sbjctLength,
                                writer.blastSpec.params,
                                writer.blastSpec.score);
    if (ret)
    {
        ::std::cerr << "Error calculating BLAST e-value and bits-score\n";
        return;
    }
    //DEBUG
//     ::std::cerr << "CONTIG index: " << alignedRead.contigId
//                 << "\t #names: " << length(writer._store.contigNameStore)
//                 << "\t #seqs: " << length(writer._store.contigStore)
//                 << '\n';
//     ::std::cerr << "READ index: " << alignedRead.readId
//                 << "\t #names: " << length(writer._store.readNameStore)
//                 << "\t #seqs: " << length(writer._store.readSeqStore)
//                 << '\n';


    const CharString sId = writer._store.contigNameStore[alignedRead.contigId];
//     const CharString qId = "foobar";
    const CharString qId = writer._store.readNameStore[alignedRead.readId];
    // Write record.[+1 added because blast is 1-index, not 0-indexed]
    writeRecord(writer._stream, qId, sId, percentIdent, ali_length, mismatches,
                gap_openings, 1, qryLength,
                alignedRead.beginPos+1, alignedRead.endPos+1, eval, bitScore,
                BlastTab(),
                typename TSpec::TProgramType());
}

// // Single-End
// template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
//           typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
// inline void onMatch(Writer<TGenome, TReads, BlastTabHdr, TDistance, TSpec> & writer,
//                     TContigId contigId,
//                     TContigPos beginPos,
//                     TContigPos endPos,
//                     TReadId readId,
//                     TErrors errors,
//                     bool reverseComplemented)
// {
//     return onMatch(writer,
//                    contigId,
//                    beginPos,
//                     TContigPos endPos,
//                     TReadId readId,
//                     TErrors errors,
//                     bool reverseComplemented)
// }
// Paired-End
template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Writer<TGenome, TReads, BlastTab, TDistance, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPosFwd,
                    TContigPos endPosFwd,
                    TReadId readIdFwd,
                    TErrors errorsFwd,
                    TContigPos beginPosRev,
                    TContigPos endPosRev,
                    TReadId readIdRev,
                    TErrors errorsRev)
{
    return;
}

// ----------------------------------------------------------------------------
// Function _alignRead()                              [Writer<HammingDistance>]
// ----------------------------------------------------------------------------
/*
//TODO(h4nn3s): implement

template <typename TGenome, typename TReads, typename TFormat, typename TSpec, typename TAlign>
inline void _alignRead(Writer<TGenome, TReads, TFormat, HammingDistance, TSpec> & writer,
                       TAlign & align,
                       TAlignedReadStoreElement & alignedRead,
                       TAlignQualityStoreElement &,
                       bool reverseComplemented)
{
    typedef TFragmentStore::TReadSeq    TReadSeq;

    resize(rows(align), 2);

    assignSource(row(align, 0), infix(writer._store.contigStore[alignedRead.contigId].seq,
                                      std::min(alignedRead.beginPos, alignedRead.endPos),
                                      std::max(alignedRead.beginPos, alignedRead.endPos)));

    TReadSeqStoreSize readId = alignedRead.readId;
    if (reverseComplemented)
        readId += (writer.reads)->readsCount;

    TReadSeq const & readSeq = writer._store.readSeqStore[readId];
    assignSource(row(align, 1), readSeq);
}

// ----------------------------------------------------------------------------
// Function _alignRead()                                 [Writer<EditDistance>]
// ----------------------------------------------------------------------------

//TODO(h4nn3s): implement

template <typename TGenome, typename TReads, typename TFormat, typename TSpec, typename TAlign>
inline void _alignRead(Writer<TGenome, TReads, TFormat, EditDistance, TSpec> & writer,
                       TAlign & align,
                       TAlignedReadStoreElement & alignedRead,
                       TAlignQualityStoreElement & alignQuality,
                       bool reverseComplemented)
{
    typedef TFragmentStore::TReadSeq    TReadSeq;

    resize(rows(align), 2);

    assignSource(row(align, 0), infix(writer._store.contigStore[alignedRead.contigId].seq,
                                      std::min(alignedRead.beginPos, alignedRead.endPos),
                                      std::max(alignedRead.beginPos, alignedRead.endPos)));

    TReadSeqStoreSize readId = alignedRead.readId;
    if (reverseComplemented)
        readId += (writer.reads)->readsCount;

    TReadSeq const & readSeq = writer._store.readSeqStore[readId];
    assignSource(row(align, 1), readSeq);

    // In this case no indels are possible.
    if ((alignQuality.errors <= 1) && (length(row(align, 0)) == length(row(align, 1))))
        return;

    globalAlignment(align, Score<short, EditDistance>(),
                    (short)-alignQuality.errors, (short)alignQuality.errors,
                    NeedlemanWunsch());
}*/

// ----------------------------------------------------------------------------
// Function _alignRead()                                 [Writer<TScore>]
// ----------------------------------------------------------------------------

//TODO(h4nn3s): implement Arbitrary scoring
/*
template <typename TGenome, typename TReads, typename TScore, typename TSpec, typename TAlign>
inline void _alignRead(Writer<TGenome, TReads, BlastTab, TScore, TSpec> & writer,
                       TAlign & align,
                       TAlignedReadStoreElement & alignedRead,
                       bool reverseComplemented)
{
    typedef TFragmentStore::TReadSeq    TReadSeq;

    resize(rows(align), 2);

    assignSource(row(align, 0), infix(writer._store.contigStore[alignedRead.contigId].seq,
                                      std::min(alignedRead.beginPos, alignedRead.endPos),
                                      std::max(alignedRead.beginPos, alignedRead.endPos)));

    TReadSeqStoreSize readId = alignedRead.readId;
//     if (reverseComplemented)
//         readId += (writer.reads)->readsCount;

    TReadSeq const & readSeq = writer._store.readSeqStore[readId];
    assignSource(row(align, 1), readSeq);


//     globalAlignment(align, writer.score, Gotoh());

    globalAlignment(align, Score<short, EditDistance>(),
                    (short)-alignQuality.errors, (short)alignQuality.errors,
                    NeedlemanWunsch());
}*/



#endif  // #ifndef SEQAN_EXTRAS_MASAI_WRITER_H_
