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
// writer.hpp: code to write out matches
// ==========================================================================

#ifndef SEQAN_SCARLET_WRITER_H_
#define SEQAN_SCARLET_WRITER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/blast.h>
#include <seqan/align.h>

#include <type_traits>

#include "options.hpp"
#include "finder.hpp"
#include "trans.hpp"
#include "misc.hpp"




template <typename TString,
          typename TScore,
          typename TStream,
          typename TSpec,
          BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
alignMatchAndWrite(TStream & stream,
                   Match const & mc,
                   StringSet<TString, TSpec> const & qrySeqs,
                   StringSet<TString, TSpec> const & subjSeqs,
                   StringSet<CharString, Owner<ConcatDirect<>>> const & qryIds,
                   StringSet<CharString, Owner<ConcatDirect<>>> const & subjIds,
                   ScarletOptions const & options,
                   TScore const & scoreScheme,
                   typename BlastStatistics<TScore>::TParams const & blastParams,
                   BlastFormat<mf, p, g> const & /*tag*/)
{
    typedef BlastFormat<mf,p,g> TFormat;

    Match m(mc); // copy match

//     // THIS SHOULD NOT BE NECESSARY; FIND OUT WHY IT IS
//     typedef  decltype(m.qryEnd) T;
//     T  qryEnd = std::min(T(length(value(qrySeqs, m.qryId)) - 1), m.qryEnd);
//     T subjEnd = std::min(T(length(value(subjSeqs, m.subjId)) - 1), m.subjEnd);


    if ((m.qryStart >= m.qryEnd) || (m.subjStart >= m.subjEnd))
    {
        std::cout << "## EARLY FAIL\n";
        _printMatch(m);

        std::cout << "QryEnd: " << m.qryEnd << "\tsubjEnd: " << m.subjEnd << std::endl;
        std::cout << "length(qry): " << length(value(qrySeqs, m.qryId))
                << "\tlength(subj): " << length(value(subjSeqs, m.subjId))
                << std::endl << std::endl;

        return 55555;
    }

//     auto qryInfix = infix(value(qrySeqs, m.qryId),
//                                 m.qryStart,
//                                 m.qryEnd);
//     auto subjInfix  = infix(value(subjSeqs, m.subjId),
//                                 m.subjStart,
//                                 m.subjEnd);

//     std::cout << "Query Id: " << m.qryId
//               << "\t TrueQryId: " << getTrueQryId(m, options, TFormat())
//               << "\t length(qryIds): " << length(qryIds)
//               << "Subj Id: " << m.subjId
//               << "\t TrueSubjId: " << getTrueSubjId(m, options, TFormat())
//               << "\t length(subjIds): " << length(subjIds) << "\n\n";

    CharString const & qryId = qryIds[getTrueQryId(m, options, TFormat())];
    CharString const & subjId = subjIds[getTrueSubjId(m, options, TFormat())];

    Align<TString, ArrayGaps> align;
    resize(rows(align), 2);
    auto & row0  = row(align, 0);
    auto & row1  = row(align, 1);

    assignSource(    row0, infix(value(qrySeqs, m.qryId),
                                 m.qryStart,
                                 m.qryEnd));
//     setBeginPosition(row0, m.qryStart);
//     setEndPosition(  row0, m.qryEnd);

    assignSource(    row1, infix(value(subjSeqs, m.subjId),
                                 m.subjStart,
                                 m.subjEnd));
//     setBeginPosition(row1, m.subjStart);
//     setEndPosition(  row1, m.subjEnd);

    std::cout << "== Positions\n";
    std::cout << "   " <<  m.qryStart << " - " << m.qryEnd << " [before ali]\n";
    int scr = 0;

    // we know the bounds
    unsigned short maxDist =      (m.qryEnd - m.qryStart)
                                / options.seedLength
                                * options.maxSeedDist;

    // compute alignment TODO extra case for hammingdistance
    if (options.semiGlobal)
    {
        // free end-gaps on subject
        AlignConfig<false, true, true, false> alignConfig;

        scr = globalAlignment(align,
                                scoreScheme,
                                alignConfig,
                                -maxDist,
                                maxDist,
                                Gotoh());
        // TODO get new positions on subj-sequence
    }
    else // local
    {
        scr = localAlignment(align,
                             scoreScheme,
                             -maxDist,
                             maxDist);
        // save new bounds of alignment
        m.qryEnd    =  m.qryStart  + endPosition(row0) - 1;
        m.qryStart  += beginPosition(row0);

        m.subjEnd   =  m.subjStart + endPosition(row1) - 1;
        m.subjStart += beginPosition(row1);

    }


    std::cout << "   " <<  m.qryStart << " - " << m.qryEnd << " [after ali]\n";
    // redefine infixes
//     qryInfix = infix(value(qrySeqs, m.qryId),
//                             m.qryStart,
//                             m.qryEnd);
//     subjInfix= infix(value(subjSeqs, m.subjId),
//                             m.subjStart,
//                             m.subjEnd);

    // extension (ungapped XDrop for now, TODO use)
    if ((!options.semiGlobal) && scr > 0)
    {
        auto max_score = scr;
        unsigned num_neg_score = 0;
        // extend left only if localAlignment didn't already shorten the match
        if (m.qryStart == mc.qryStart)
        {

            while ((m.qryStart > 0) &&
                   (m.subjStart > 0) &&
                   ((max_score - scr) >= options.xDropOff))
            {
                --m.qryStart;
                --m.subjStart;
                scr += score(scoreScheme,
                               sequenceEntryForScore(scoreScheme,
                                                     row0,
                                                     m.qryStart),
                               sequenceEntryForScore(scoreScheme,
                                                     row1,
                                                     m.subjStart));
                if (scr >= max_score)
                {
                    max_score = scr;
                    num_neg_score = 0;
                }
                else
                {
                    ++num_neg_score;
                }
            }
            // reset to max's position
            scr = max_score;
            m.qryStart  += num_neg_score;
            m.subjStart += num_neg_score;
        }

        max_score = scr;
        num_neg_score = 0;
        // extend right only if localAlignment didn't already shorten the match
        if (m.qryEnd == mc.qryEnd)
        {

            while ((m.qryEnd > 0) &&
                   (m.subjEnd > 0) &&
                   ((max_score - scr) >= options.xDropOff))
            {
                ++m.qryEnd;
                ++m.subjEnd;
                scr += score(scoreScheme,
                               sequenceEntryForScore(scoreScheme,
                                                     row0,
                                                     m.qryEnd),
                               sequenceEntryForScore(scoreScheme,
                                                     row1,
                                                     m.subjEnd));
                if (scr >= max_score)
                {
                    max_score = scr;
                    num_neg_score = 0;
                }
                else
                {
                    ++num_neg_score;
                }
            }
            // reset to max's position
            scr = max_score;
            m.qryEnd  -= num_neg_score;
            m.subjEnd -= num_neg_score;
        }

        typename Value<decltype(row0._source)>::Type newseq0;
        if (mc.qryStart - m.qryStart != 0)
        {
            row0._array[0] += mc.qryStart - m.qryStart;
            append(newseq0, infix(value(qrySeqs, m.qryId),
                                    mc.qryStart,
                                    m.qryStart));
        }
        append(newseq0, value(row0._source));
        if (m.qryEnd - mc.qryEnd != 0)
        {
            row0._array[length(row0._array)-1] += m.qryEnd - mc.qryEnd;
            append(newseq0, infix(value(qrySeqs, m.qryId),
                                    m.qryEnd,
                                    mc.qryEnd));
        }
        setValue(row0._source, newseq0);


        typename Value<decltype(row1._source)>::Type newseq1;
        if (mc.subjStart - m.subjStart > 0 )
        {
            if (isGap(row1, 0))
                insertValue(row1._array, 0, mc.subjStart - m.subjStart);
            else
                row1._array[0] += mc.subjStart - m.subjStart;
            append(newseq1, infix(value(subjSeqs, m.subjId),
                                    mc.subjStart,
                                    m.subjStart));
        }
        append(newseq1, value(row1._source));
        if (mc.subjEnd - m.subjEnd)
        {
            if (isGap(row1, length(row1)-1))
                appendValue(row1._array, m.subjEnd - mc.subjEnd);
            else
                row1._array[length(row1._array)-1] += m.subjEnd - mc.subjEnd;

            append(newseq1, infix(value(subjSeqs, m.subjId),
                                    m.subjEnd,
                                    mc.subjEnd));
        }
        setValue(row1._source, newseq1);

        std::cout << "   " <<  m.qryStart << " - " << m.qryEnd << " [after ext]\n";
    }

    if (scr <= 0)
    {
        std::cout << "## LATE FAIL\n";
        _printMatch(mc);
         CharString q = infix(value(qrySeqs, mc.qryId),
                                 mc.qryStart,
                                 mc.qryEnd);
         CharString s = infix(value(subjSeqs, mc.subjId),
                                 mc.subjStart,
                                 mc.subjEnd);
        std::cout << "QRY:  " << toCString(q) << '\n'
                  << "SUBJ: " << toCString(s) << "\n\n";
        _printMatch(m);


        std::cout << "length(qry): " << length(value(qrySeqs, m.qryId))
                << "\tlength(subj): " << length(value(subjSeqs, m.subjId))
                << std::endl;

        return 54321;
    }
//     std::cout << "##LINE: " << __LINE__ << '\n';
    long sc = 0;
    unsigned int ali_length = 0;
    unsigned int identities = 0;
    unsigned int positives = 0;
    unsigned int mismatches = 0;
    unsigned int gaps = 0;
    unsigned int gap_openings = 0;

    int ret = calcStatsAndScore(sc, ali_length, identities,
                                positives,mismatches, gaps,
                                gap_openings, row0, row1,
                                scoreScheme);
    if (ret)
    {
        ::std::cerr << "Error calculating BLAST stats and score\n";
        return -1;
    }

    double eval = 0;
    double bitScore = 0;
    const  unsigned long sbjctLength = length(value(subjSeqs, m.subjId));
    const unsigned long qryLength = length(value(qrySeqs, m.qryId));
    ret = calcBitScoreAndEValue(bitScore, eval,
                                sc, sbjctLength, qryLength,
                                blastParams,
                                scoreScheme);
    if (ret)
    {
        ::std::cerr << "Error calculating BLAST e-value and bits-score\n";
        return -1;
    }

    if (eval > options.eCutOff)
    {
//         CharString s1 = qryInfix;
//         CharString s2 = subjInfix;
//         std::cout << "## Hit discarded with:\n"
//                   << "  qry: " << toCString(s1)
//                   << "\n subj: " << toCString(s2)
//                   << "\n scores::\tGotoh " << score
//                   << "\tblast " << sc
//                   << "\tbits " << bitScore
//                   << "\tevalue " << eval << "\n";
        return 12345;
    }
    // blast is 1-index, not 0-indexed
    unsigned realQryBeginPos  = getTrueQryStartPos (m, options, TFormat()) + 1;
    unsigned realQryEndPos    = getTrueQryEndPos   (m, options, TFormat()) + 1;
    unsigned realSubjBeginPos = getTrueSubjStartPos(m, options, TFormat()) + 1;
    unsigned realSubjEndPos   = getTrueSubjEndPos  (m, options, TFormat()) + 1;

    // TODO lockguard this so the loop around this can be parallized
    ret = writeRecord(stream,
                      qryId, subjId,
                      identities, ali_length, mismatches, gap_openings,
                      realQryBeginPos, realQryEndPos,
                      realSubjBeginPos, realSubjEndPos,
                      eval, bitScore, TFormat());
    return ret;
}




template <typename TString,
          typename TSpec,
          typename TScore,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
iterateMatches(StringSet<TString, TSpec> const & qrySeqs,
             StringSet<TString, TSpec> const & subjSeqs,
             String<Match> const & ms,
             ScarletOptions const & options,
             TScore const & scoreScheme,
             typename BlastStatistics<TScore>::TParams const & blastParams,
             BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    std::ofstream stream;
    stream.open(toCString(options.output));

    double start = sysTime();
    std::cout << "Loading Qry Ids..." << std::flush;

    StringSet<CharString, Owner<ConcatDirect<> > > qryIds;
    int ret = 0;
    if (options.fileFormat == 0)
        ret = loadIds(qryIds, options.queryFile, Fasta());
    else
        ret = loadIds(qryIds, options.queryFile, Fastq());
    if (ret)
        return ret;

    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

    start = sysTime();
    std::cout << "Loading Subj Ids..." << std::flush;

    StringSet<CharString, Owner<ConcatDirect<> > > subjIds;
    ret = loadIds(subjIds, options.dbFile, Fasta());
    if (ret)
        return ret;

    std::cout << " done.\n";
    finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

    unsigned long goodMatches = 0; // eval <= 0.1
    unsigned long badMatches = 0;  // eval > 0.1
    unsigned long earlyNotMatches = 0;  // something went wrong this aint a match
    unsigned long lateNotMatches = 0;  // something went wrong this aint a match
    for (Match const & ma : ms)
    {
        ret = alignMatchAndWrite(stream, ma, qrySeqs, subjSeqs, qryIds, subjIds,
                                 options, scoreScheme, blastParams, TFormat());
        if (ret == 0)
            ++goodMatches;
        else if (ret == 12345)
            ++badMatches;
        else if (ret == 54321)
            ++lateNotMatches;
        else if (ret == 55555)
            ++earlyNotMatches;
        else
            return ret;
    }
    std::cout << "NUMBER OF MATCHES:\n" << "Good: " << goodMatches
    << "\nBad: " << badMatches
    << "\nFailed early: " << earlyNotMatches << "\n"
    << "\nFailed late: " << lateNotMatches << "\n";
    return 0;
}

template <typename TString,
          typename TSpec,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline int
prepareScoring(StringSet<TString, TSpec> const & qrySeqs,
               StringSet<TString, TSpec> const & dbSeqs,
               String<Match>             const & ms,
               ScarletOptions            const & options,
               BlastFormat<m,
                           BlastFormatOptions::BlastN,
                           g>            const & /*tag*/)
{
    typedef BlastFormat<m,BlastFormatOptions::BlastN,g> TFormat;

    Score<int, Simple> scoreScheme(2, -3, -2, -5); // default blastn

    typename BlastStatistics<Score<int, Simple> >::TParams blastParams;

    int ret = getScoringParams(blastParams, scoreScheme);
    if (ret)
    {
        ::std::cerr << "Could not computer Karlin-Altschul-Values for "
                    << "Scoring Scheme. Exiting.\n";
        return -1;
    }

    return iterateMatches(qrySeqs, dbSeqs, ms, options,
                            scoreScheme, blastParams, TFormat());
}

template <typename TString,
          typename TSpec,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
prepareScoring(StringSet<TString, TSpec> const & qrySeqs,
               StringSet<TString, TSpec> const & dbSeqs,
               String<Match>             const & ms,
               ScarletOptions            const & options,
               BlastFormat<m, p, g>      const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    Blosum62    scoreScheme;
    setScoreGapOpen(scoreScheme, -11);
    setScoreGapExtend(scoreScheme, -1);

    typename BlastStatistics<Blosum62>::TParams blastParams;

    int ret = getScoringParams(blastParams, scoreScheme);
    if (ret)
    {
        ::std::cerr << "Could not computer Karlin-Altschul-Values for "
                    << "Scoring Scheme. Exiting.\n";
        return -1;
    }
    return iterateMatches(qrySeqs, dbSeqs, ms, options,
                          scoreScheme, blastParams, TFormat());
}



// load the query sequences in (possibly) translated, but
// unreduced alphabet

template <typename TString,
          typename TSpec,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
loadQrySequences(StringSet<TString, TSpec> const & dbSeqs,
                 String<Match>             const & ms,
                 ScarletOptions            const & options,
                 BlastFormat<m, p, g>      const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    typedef typename UnreducedStringSet<p>::Type TQrySeqs;

    double start = sysTime();
    std::cout << "Loading Qry Sequences..." << std::flush;

    TQrySeqs qrySeqs;
    CharString _qrySeqs = options.queryFile;
    append(_qrySeqs, ".unredqry");
    int ret = open(qrySeqs, toCString(_qrySeqs));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }

    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "Amount: " << length(qrySeqs) << "\n\n"<< std::flush;

    return prepareScoring(qrySeqs, dbSeqs, ms, options, TFormat());
}

// load the database sequences in (possibly) translated, but
// unreduced alphabet
template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
loadDbSequences(String<Match>        const & ms,
                ScarletOptions       const & options,
                BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    double start = sysTime();
    std::cout << "Loading Subj Sequences..." << std::flush;

    typedef typename UnreducedStringSet<p>::Type TSubjSeqs;
    TSubjSeqs dbSeqs;
    CharString _dbSeqs = options.dbFile;
    if (options.alphReduction > 0)
        append(_dbSeqs, ".unredsubj"); // get unreduced stringset
    else // stringset already dumped by index dump
        append(_dbSeqs, ".txt");
    int ret = open(dbSeqs, toCString(_dbSeqs));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }

    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "Amount: " << length(dbSeqs) << "\n\n"<< std::flush;

    return loadQrySequences(dbSeqs, ms, options, TFormat());
}



template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writePipeline(String<Match>        const & ms,
              ScarletOptions       const & options,
              BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;
    return loadDbSequences(ms, options, TFormat());
}

#endif // header guard