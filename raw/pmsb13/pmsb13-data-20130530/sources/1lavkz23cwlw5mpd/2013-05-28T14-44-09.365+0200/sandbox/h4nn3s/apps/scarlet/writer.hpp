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
#include <seqan/seeds.h>

#include <type_traits>

#include "options.hpp"
#include "finder.hpp"
#include "trans.hpp"
#include "misc.hpp"




template <typename TString,
          typename TScore,
          typename TSpec,
          BlastFormatOptions::M mf,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
computeBlastMatch(BlastMatch & bm,
                   StringSet<TString, TSpec> const & qrySeqs,
                   StringSet<TString, TSpec> const & subjSeqs,
//                    StringSet<CharString, Owner<ConcatDirect<>>> const & qryIds,
//                    StringSet<CharString, Owner<ConcatDirect<>>> const & subjIds,
                   ScarletOptions const & options,
                   TScore const & scoreScheme,
                   typename BlastStatistics<TScore>::TParams const & blastParams,
                   BlastFormat<mf, p, g> const & /*tag*/)
{
    typedef BlastFormat<mf,p,g> TFormat;


//     // THIS SHOULD NOT BE NECESSARY; FIND OUT WHY IT IS
//     typedef  decltype(m.qryEnd) T;
//     T  qryEnd = std::min(T(length(value(qrySeqs, bm.m.qryId)) - 1), bm.m.qryEnd);
//     T subjEnd = std::min(T(length(value(subjSeqs, bm.m.subjId)) - 1), bm.m.subjEnd);


    if ((bm.m.qryStart >= bm.m.qryEnd) || (bm.m.subjStart >= bm.m.subjEnd))
    {
        std::cout << "## EARLY FAIL\n";
        _printMatch(bm.m);

        std::cout << "QryEnd: " << bm.m.qryEnd << "\tsubjEnd: " << bm.m.subjEnd << std::endl;
        std::cout << "length(qry): " << length(value(qrySeqs, bm.m.qryId))
                << "\tlength(subj): " << length(value(subjSeqs, bm.m.subjId))
                << std::endl << std::endl;

        return 55555;
    }

//     auto qryInfix = infix(value(qrySeqs, bm.m.qryId),
//                                 bm.m.qryStart,
//                                 bm.m.qryEnd);
//     auto subjInfix  = infix(value(subjSeqs, bm.m.subjId),
//                                 bm.m.subjStart,
//                                 bm.m.subjEnd);

//     std::cout << "Query Id: " << bm.m.qryId
//               << "\t TrueQryId: " << getTrueQryId(bm.m, options, TFormat())
//               << "\t length(qryIds): " << length(qryIds)
//               << "Subj Id: " << bm.m.subjId
//               << "\t TrueSubjId: " << getTrueSubjId(bm.m, options, TFormat())
//               << "\t length(subjIds): " << length(subjIds) << "\n\n";

    Align<TString, ArrayGaps> align;
    resize(rows(align), 2);
    auto & row0  = row(align, 0);
    auto & row1  = row(align, 1);

    assignSource(row0, infix(value(qrySeqs, bm.m.qryId),
                             bm.m.qryStart,
                             bm.m.qryEnd));
    assignSource(row1, infix(value(subjSeqs, bm.m.subjId),
                             bm.m.subjStart,
                             bm.m.subjEnd));

    std::cout << "== Positions\n";
    std::cout << "   " <<  bm.m.qryStart << " - " << bm.m.qryEnd << " [before ali]\n";
    std::cout << align << std::endl;
    int scr = 0;

    // we know the bounds
    unsigned short maxDist =    (bm.m.qryEnd - bm.m.qryStart)
                              / (options.seedLength - 1)
                              *  options.maxSeedDist;

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
        bm.m.subjEnd   =   bm.m.subjStart
                      + toSourcePosition(row1,
                                         toViewPosition(row0, length(row0) -1));
        bm.m.subjStart +=  toSourcePosition(row1,
                                         toViewPosition(row0, 0));
    }
    else // local
    {
        scr = localAlignment(align,
                             scoreScheme,
                             -maxDist,
                             maxDist);
        // save new bounds of alignment
        bm.m.qryEnd    =  bm.m.qryStart  + endPosition(row0) - 1;
        bm.m.qryStart  += beginPosition(row0);

        bm.m.subjEnd   =  bm.m.subjStart + endPosition(row1) - 1;
        bm.m.subjStart += beginPosition(row1);

    }


    std::cout << "   " <<  bm.m.qryStart << " - " << bm.m.qryEnd << " [after ali]\n";
    std::cout << align << std::endl;

    if ((!options.semiGlobal) && (scr > 0))
    {

        Seed<Simple>           seed(bm.m.subjStart, bm.m.qryStart,
                                    bm.m.subjEnd, bm.m.qryEnd);
        extendSeed(seed,
                   value(subjSeqs, bm.m.subjId),
                   value(qrySeqs, bm.m.qryId),
                   EXTEND_BOTH,
                   scoreScheme,
                   int(options.xDropOff),
                   UnGappedXDrop());
        bm.m.subjStart = beginPositionH(seed);
        bm.m.qryStart  = beginPositionV(seed);
        bm.m.subjEnd   = endPositionH(seed);
        bm.m.qryEnd    = endPositionV(seed);

        // TODO: instead of global realigned, align only new end and beginning
        // and concatenate with old aling object
        assignSource(row0, infix(value(qrySeqs, bm.m.qryId),
                                bm.m.qryStart,
                                bm.m.qryEnd));
        assignSource(row1, infix(value(subjSeqs, bm.m.subjId),
                                bm.m.subjStart,
                                bm.m.subjEnd));
        // true global Alignment
        AlignConfig<false, false, false, false> alignConfig;

        scr = globalAlignment(align,
                              scoreScheme,
                              alignConfig,
                              Gotoh());

        std::cout << "   " <<  bm.m.qryStart << " - " << bm.m.qryEnd << " [after ext]\n";
        std::cout << align << std::endl;
    }

    if (scr <= 0)
    {
        std::cout << "## LATE FAIL\n";
        _printMatch(bm.orig_m);
         CharString q = infix(value(qrySeqs, bm.orig_m.qryId),
                                 bm.orig_m.qryStart,
                                 bm.orig_m.qryEnd);
         CharString s = infix(value(subjSeqs, bm.orig_m.subjId),
                                 bm.orig_m.subjStart,
                                 bm.orig_m.subjEnd);
        std::cout << "QRY:  " << toCString(q) << '\n'
                  << "SUBJ: " << toCString(s) << "\n\n";
        _printMatch(bm.m);


        std::cout << "length(qry): " << length(value(qrySeqs, bm.m.qryId))
                << "\tlength(subj): " << length(value(subjSeqs, bm.m.subjId))
                << std::endl;

        return 54321;
    }
//     std::cout << "##LINE: " << __LINE__ << '\n';
    long sc = 0;

    int ret = calcStatsAndScore(sc, bm.ali_length, bm.identities,
                                bm.positives,bm.mismatches, bm.gaps,
                                bm.gap_openings, row0, row1,
                                scoreScheme);
    if (ret)
    {
        ::std::cerr << "Error calculating BLAST stats and score\n";
        return -1;
    }

    const  unsigned long sbjctLength = length(value(subjSeqs, bm.m.subjId));
    const unsigned long qryLength = length(value(qrySeqs, bm.m.qryId));
    ret = calcBitScoreAndEValue(bm.bitScore, bm.eval,
                                sc, sbjctLength, qryLength,
                                blastParams,
                                scoreScheme);
    if (ret)
    {
        ::std::cerr << "Error calculating BLAST e-value and bits-score\n";
        return -1;
    }

    if (bm.eval > options.eCutOff)
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

    // UNTRANSLATE
    // blast is 1-index, not 0-indexed
    bm.m.qryId     = getTrueQryId       (bm.m, options, TFormat());
    bm.m.qryStart  = getTrueQryStartPos (bm.m, options, TFormat()) + 1;
    bm.m.qryEnd    = getTrueQryEndPos   (bm.m, options, TFormat()) + 1;
    bm.m.subjId    = getTrueSubjId      (bm.m, options, TFormat());
    bm.m.subjStart = getTrueSubjStartPos(bm.m, options, TFormat()) + 1;
    bm.m.subjEnd   = getTrueSubjEndPos  (bm.m, options, TFormat()) + 1;

    // TODO if TFormat == Report then add align to bm

    return 0;
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


    // partition matches into intervals per query sequence
    // this only works since matches are sorted by qryId AND
    // different qryIds that map to the same TrueQryId appear consecutively
    unsigned long qrysWithMatches = 0;
    for (unsigned long l = 0; l < length(ms); ++l)
    {
        if ((l>0) &&
            (getTrueQryId(ms[l-1], options, TFormat())
             != getTrueQryId(ms[l], options, TFormat())))
            ++qrysWithMatches;
    }
    std::vector<std::tuple<Match::TQId,   // TrueQryId
                           unsigned long, // begin of match interval
                           unsigned long> // end of match interval
                           > intervals;
    intervals.resize(qrysWithMatches);
    for (unsigned long l = 0, l2 = 0; l < length(ms); ++l)
    {
        if (l == 0)
        {
            std::get<0>(intervals[l2]) = getTrueQryId(ms[l], options, TFormat());
            std::get<1>(intervals[l2]) = l;
        }
        // new query
        else if (std::get<0>(intervals[l2]) != getTrueQryId(ms[l], options, TFormat()))
        {
            // set end position of last
            std::get<2>(intervals[l2]) = l-1; 

            ++l2;
            std::get<0>(intervals[l2]) = getTrueQryId(ms[l], options, TFormat());
            std::get<1>(intervals[l2]) = l;
        }
    }
    // TODO try guided instead
    int error = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (unsigned long l = 0; l < length(intervals); ++l)
    {
        #pragma omp flush(error)
        if (error == 0)
        {
            auto const & in = intervals[l];

            std::list<BlastMatch> bmatches;
            //TODO inner loop not parallelized right now
            for (unsigned long beg = std::get<1>(in);
                beg < std::get<2>(in);
                ++beg)
            {
                BlastMatch bm(ms[l]);
                int ret = computeBlastMatch(bm, qrySeqs, subjSeqs, options,
                                            scoreScheme, blastParams, TFormat());

                if (ret == 0)
                    bmatches.push_back(bm);

                if (ret == 0)
                {
                    #pragma omp atomic
                    ++goodMatches;
                } else if (ret == 12345)
                {
                    #pragma omp atomic
                    ++badMatches;
                } else if (ret == 54321)
                {
                    #pragma omp atomic
                    ++lateNotMatches;
                } else if (ret == 55555)
                {
                    #pragma omp atomic
                    ++earlyNotMatches;
                } else
                {
                    #pragma omp critical
                    error = ret;
                }
            }

            bmatches.sort();
            //TODO remove duplicates

            #pragma omp critical
            {
                for (BlastMatch const & bm : bmatches)
                {
                    int ret = writeRecord(stream,
                                        qryIds[bm.m.qryId],
                                        subjIds[bm.m.subjId],
                                        bm.identities,
                                        bm.ali_length,
                                        bm.mismatches,
                                        bm.gap_openings,
                                        bm.m.qryStart,
                                        bm.m.qryEnd,
                                        bm.m.subjStart,
                                        bm.m.subjEnd,
                                        bm.eval,
                                        bm.bitScore,
                                        TFormat());
                    if (ret)
                    {
                        error = ret;
                        break;
                    }
                }
            }
        }
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