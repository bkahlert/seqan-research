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
// scarlet_indexer.hpp: Main File for the indexer application
// ==========================================================================

#ifndef SEQAN_SCARLET_SCARLET_INDEXER_H_
#define SEQAN_SCARLET_SCARLET_INDEXER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/index_extras.h>

#include "misc.hpp"
#include "options.hpp"
#include "trans.hpp"
#include "alph.hpp"


using namespace seqan;

template <typename TString, typename TSpec,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step03_generateIndexAndDump(StringSet<TString, TSpec> & seqs,
                            ScarletIndexerOptions const & options,
                            BlastFormat<m, p, g> const & /*tag*/)
{
    typedef Index<StringSet<TString, TSpec>, IndexSa<> > TDbIndex;

    // Generate Index
    std::cout << "Generating Index..." << std::flush;
    double s = sysTime();
    TDbIndex dbIndex(seqs);
    typename Iterator<TDbIndex, TopDown<> >::Type it(dbIndex); // instantiate
    double e = sysTime() - s;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << e << "s \n\n" << std::flush;


    // Dump Index
    std::cout << "Writing Index to disk..." << std::flush;
    s = sysTime();
    save(dbIndex, toCString(options.dbFile));
    e = sysTime() - s;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << e << "s \n" << std::flush;

    return 0;
}

template <typename TString, typename TSpec,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step02_reduceAlphabet(StringSet<TString, TSpec> & oldSubj,
                      ScarletIndexerOptions const & options,
                      BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    if (options.alphReduction == 10)
    {
        StringSet<String<AminoAcid10>, Owner<ConcatDirect<>>> newSubj;
        newSubj.concat = oldSubj.concat; //implicit conversion
        newSubj.limits = oldSubj.limits;
        return step03_generateIndexAndDump(newSubj, options, TFormat());
    }

    return step03_generateIndexAndDump(oldSubj, options, TFormat());
}

// --------------------------------------------------------------------------
// Function preprocessSubjSeqs()
// --------------------------------------------------------------------------




template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                   TOldSubj & oldSubj,
                                   ScarletIndexerOptions const & /**/,
                                   BlastFormat<m,
                                               BlastFormatOptions::BlastN,
                                               g> const & /*tag*/)
{
    newSubj.concat = oldSubj.concat; //implicit conversion
    newSubj.limits = oldSubj.limits;
}

template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                    TOldSubj const & oldSubj,
                                    ScarletIndexerOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::BlastP,
                                                g> const & /*tag*/)
{
    newSubj.concat = oldSubj.concat; //implicit conversion
    newSubj.limits = oldSubj.limits;

}

template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                    TOldSubj const & oldSubj,
                                    ScarletIndexerOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::BlastX,
                                                g> const & /*tag*/)
{
    newSubj.concat = oldSubj.concat; //implicit conversion
    newSubj.limits = oldSubj.limits;
}

template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                    TOldSubj const & oldSubj,
                                    ScarletIndexerOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::TBlastN,
                                                g> const & /*tag*/)
{
    reserve(newSubj.concat, length(oldSubj.concat) / 3 * 6, Exact());
    //  / 3 * 6 is a little too much, but the speedup is worth it
    for (auto const & s : oldSubj)
        appendSixFrameTranslation(newSubj, s);
}

template <typename TNewSubj,
          typename TOldSubj,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step01a_preprocessSubj(TNewSubj & newSubj,
                                    TOldSubj const & oldSubj,
                                    ScarletIndexerOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::TBlastX,
                                                g> const & /*tag*/)
{
    reserve(newSubj.concat, length(oldSubj.concat) / 3 * 6, Exact());
    //  / 3 * 6 is a little too much, but the speedup is worth it
    for (auto const & s : oldSubj)
        appendSixFrameTranslation(newSubj, s);
}

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step01_preprocessSubjSeqs(
                       StringSet<CharString, Owner<ConcatDirect<> > >  & oldSubj,
                       ScarletIndexerOptions       const & options,
                       BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    typename UnreducedStringSet<p>::Type newSubj;

    double start = sysTime();
    std::cout << "Preprocessing Subj Sequences..." << std::flush;

    // depending on blastProgram
    step01a_preprocessSubj(newSubj, oldSubj, options, TFormat());

    clear(oldSubj);
    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "Number of queries after preproc: " << length(newSubj)
                << "\n\n" << std::flush;

    if (options.alphReduction > 0)
    {
        start = sysTime();
        std::cout << "Dumping unreduced Subj Sequences..." << std::flush;

        //TODO save to TMPDIR instead
        CharString _path = options.dbFile;
        append(_path, ".unredsubj");
        save(newSubj, toCString(_path));

        std::cout << " done.\n";
        finish = sysTime() - start;
        std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

        return step02_reduceAlphabet(newSubj, options, TFormat());
    }

    return step03_generateIndexAndDump(newSubj, options, TFormat());
}

// --------------------------------------------------------------------------
// Function loadSubjSequences()
// --------------------------------------------------------------------------



template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step00_loadSubjSequences(ScarletIndexerOptions const & options,
                         BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    // load sequences as CharString first to be error tolerant
    StringSet<CharString, Owner<ConcatDirect<> > > seqs;
    double start = sysTime();
    std::cout << "Loading Database Sequences..." << std::flush;
    int ret = 0;
    if (options.fileFormat == 1)
        ret = loadSequences(seqs, options.dbFile, Fastq());
    else
        ret = loadSequences(seqs, options.dbFile, Fasta());
    if (ret)
        return ret;
    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;

    unsigned long maxLen = 0ul;
    for (auto const & s : seqs)
        if (length(s) > maxLen)
            maxLen = length(s);
    std::cout << "Number of sequences read: " << length(seqs)
              << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;

    return step01_preprocessSubjSeqs(seqs, options, TFormat());
}

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
beginPipeline(ScarletIndexerOptions const & options,
              BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;
    return step00_loadSubjSequences(options, TFormat());
}


#endif // header guard