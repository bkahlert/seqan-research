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


using namespace seqan;

template <typename TString, typename TSpec>
inline int
step03_generateIndexAndDump(StringSet<TString, TSpec> & seqs,
                      ScarletIndexerOptions const & options)
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

template <typename TString, typename TSpec>
inline int
step02_reduceAlphabet(StringSet<TString, TSpec> & seqs,
                      ScarletIndexerOptions const & options)
{
    // TODO reduce alphabet
    return step03_generateIndexAndDump(seqs, options);
}


inline int
step01_preprocessDbSequences(StringSet<CharString, Owner<ConcatDirect<>>> seqs,
                             ScarletIndexerOptions const & options)
{
    double start = sysTime();
    std::cout << "Preprocessing Db Sequences..." << std::flush;

    switch (options.blastProg)
    {
        case BlastFormatOptions::BlastN : // untranslated nucleotide
        {
            StringSet<Dna5String, Owner<ConcatDirect<> > > newSeqs;
            reserve(newSeqs.concat, length(seqs.concat), Exact());
            for (auto const & s : seqs)
            {
                appendValue(newSeqs, s, Exact());
            }
            clear(seqs);
            std::cout << " done.\n";
            double finish = sysTime() - start;
            std::cout << "Runtime: " << finish << "s \n" << std::flush;
            std::cout << "Number of queries after preproc: " << length(newSeqs)
                      << "\n\n" << std::flush;
            return step03_generateIndexAndDump(newSeqs, options);
        } break;

        case BlastFormatOptions::BlastP :
        case BlastFormatOptions::BlastX : // "untranslated" protein
        {
            StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > newSeqs;
            reserve(newSeqs.concat, length(seqs.concat), Exact());
            for (auto const & s : seqs)
            {
                appendValue(newSeqs, s, Exact());
            }
            clear(seqs);
            std::cout << " done.\n";
            double finish = sysTime() - start;
            std::cout << "Runtime: " << finish << "s \n" << std::flush;
            std::cout << "Number of queries after preproc: " << length(newSeqs)
                      << "\n\n" << std::flush;
            if (options.alphReduction >= 0)
                return step02_reduceAlphabet(newSeqs, options);
            else
                return step03_generateIndexAndDump(newSeqs, options);
        } break;

        case BlastFormatOptions::TBlastN :
        case BlastFormatOptions::TBlastX : // translated nucleotide
        {
            //TODO translate
            StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > newSeqs;
            reserve(newSeqs.concat, length(seqs.concat), Exact());
            for (auto const & s : seqs)
            {
                appendValue(newSeqs, s, Exact());
            }
            clear(seqs);
            std::cout << " done.\n";
            double finish = sysTime() - start;
            std::cout << "Runtime: " << finish << "s \n" << std::flush;
            std::cout << "Number of queries after preproc: " << length(newSeqs)
                      << "\n\n" << std::flush;
            if (options.alphReduction >= 0)
                return step02_reduceAlphabet(newSeqs, options);
            else
                return step03_generateIndexAndDump(newSeqs, options);
        } break;
        default:
            return -1;
    }
    return -1;
}

inline int
step00_loadDbSequences(ScarletIndexerOptions const & options)
{
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

    return step01_preprocessDbSequences(seqs, options);
}


inline int
beginPipeline(ScarletIndexerOptions const & options)
{
    return step00_loadDbSequences(options);
}


#endif // header guard