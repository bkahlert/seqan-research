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
// scarlet.hpp: contains the main progam pipeline
// ==========================================================================


#ifndef SEQAN_SCARLET_SCARLET_H_
#define SEQAN_SCARLET_SCARLET_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index.h>
#include <seqan/index_extras.h>

#include <seqan/blast.h>

#include "options.hpp"
#include "finder.hpp"
#include "misc.hpp"
#include "trans.hpp"
#include "alph.hpp"

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function seedAndSearch()
// --------------------------------------------------------------------------

template <typename TStringSet>
struct Comp
{
//     TSeedIndex const & index;
//     typedef typename Fibre<TSeedIndex, EsaSA>::Type TSa;
//     typedef typename SAValue<TSa>::Type TSav;
    
    TStringSet const & stringSet;

    Comp(TStringSet const & _stringSet) : stringSet(_stringSet)
    {}

    inline bool operator() (unsigned long i, unsigned long j) 
    {
        return (value(stringSet,i) < value(stringSet,i));
    }
};


template <typename TString,
          typename TSpec,
          typename TDbIndex,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step04_seedAndSearch(MatchManager               & matchManager,
                     StringSet<TString, TSpec > & qrySeqs,
                     TDbIndex                   & dbIndex,
                     ScarletOptions       const & options,
                     BlastFormat<m, p, g> const & /*tag*/)
{
//     typedef BlastFormat<m,p,g> TFormat;

    // SEEDING
//     typedef StringSet<typename Infix<TString const>::Type>  TSeeds;
    typedef StringSet<CharString>  TSeeds;
    typedef Index<TSeeds, IndexSa<> >                       TSeedIndex;

    TSeeds seeds;
    std::cout << "Generating Seeds..." << std::flush;
    double start = sysTime();
    for (unsigned long i = 0; i < length(qrySeqs); ++i)
    {
//         TString const & str = value(qrySeqs, i);

        for (unsigned j = 0;
             ((j+1) * options.seedLength) - 1 < length(value(qrySeqs, i));
             ++j)
        {
//             CharString cs = infix(value(qrySeqs, i),
//                                      j* options.seedLength,
//                                      ((j+1) * options.seedLength) - 1);
//             std::cout << "I  Seed " << length(seeds) +1
//                       << '\t' << i << '\t' << j << ": "
//                       << toCString(cs) << "\n";

            appendValue(seeds, infix(value(qrySeqs, i),
                                     j* options.seedLength,
                                     ((j+1) * options.seedLength) - 1),
                        Generous());
            appendValue(matchManager.seedRefs,  i, Generous());
            appendValue(matchManager.seedRanks, j, Generous());
        }
    }
    double finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "Number of seeds created: " << length(seeds) << "\n\n";

//     unsigned j = 0;
//     for (auto const & s : seeds)
//     {
//         CharString cs = s;
//         std::cout << "II Seed " << j++ << '\t' << matchManager.seedRefs[j]
//                   << '\t' << matchManager.seedRanks[j] << ": "
//                   << toCString(cs) << "\n";
//     }

    // BUILD INDEX
    std::cout << "Generating Query-Index..." << std::flush;
    start = sysTime();
    TSeedIndex seedIndex(seeds);
    // we only want full length seed sequences in index
    {
        typedef typename Fibre<TSeedIndex, EsaSA>::Type TSa;
        TSa & sa = indexSA(seedIndex);
        
        std::vector<unsigned long> refs;
        refs.resize(length(seeds));
        for (unsigned long l = 0; l < length(refs); ++l)
            refs[l] = l;
        Comp<TSeeds> comp(seeds);
        std::sort(refs.begin(), refs.end(), comp);

        resize(sa, length(seeds));

        unsigned long count = 0;

        for (auto & sav : sa)
        {
            setSeqNo(sav, ref[count]);
            setSeqOffset(sav, 0);
            ++count;
        }
         typename Iterator<TSeedIndex, TopDown<> >::Type it(seedIndex); // instantiate
    }
    
/*
    typename Iterator<TSeedIndex, TopDown<> >::Type it(seedIndex); // instantiate
    // filter out non-full length suffices
    TSa newSa;
    for (int i = 0; i < length(sa); ++i) //TODO could this be sped up?
    {
//         std::cout << "SA\tseqNo: " << getSeqNo(value(sa,i))
//                   << "\tseqOffset: " << getSeqOffset(value(sa,i))
//                   << "\tseq; "<< suffix(seeds[getSeqNo(value(sa,i))],
//                                         getSeqOffset(value(sa,i))) << "\n";
        if (getSeqOffset(value(sa,i)) == 0)
            appendValue(newSa, value(sa,i));
    }
    sa = newSa;
    clear(newSa); */
    finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "Number of fibres in SeedIndex: " << length(sa) << "\n\n";

//     unsigned j = 0;
//     for (auto const & sav :  sa)
//     {
//         CharString cs = seeds[getSeqNo(sav)];
//         std::cout << "Fibre " << j++ << ": "
//                   << toCString(cs) << "\n";
//     }
//     for (auto & sav : indexSA(seedIndex))
//     {
//         if (prefix(dbSeqs[getSeqNo(sav)], CUTLEN)
//                 == qrySeqs[0])
//             std::cout << "juhu " << getSeqNo(sav) << '\n';
//     }

    // FIND
    std::cout << "Starting a search..." << std::flush;
    start = sysTime();
    typedef ext::Backtracking<EditDistance,
                              ext::BacktrackingSemiGlobal>  BackSpec;
    typedef ext::Finder<TDbIndex, TSeedIndex,
                        MatchManager, BackSpec>             ScarletFinder;

    ScarletFinder finder(matchManager);

    ext::find(finder, dbIndex, seedIndex, options.maxSeedDist);
    finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

    return 0;
}

// --------------------------------------------------------------------------
// Function reduceAlphabet()
// --------------------------------------------------------------------------


template <typename TDbIndex,
          typename TString, typename TSpec,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step03_reduceAlphabet(MatchManager               & matchManager,
                      StringSet<TString,TSpec>   & oldQry,
                      TDbIndex                   & dbIndex,
                      ScarletOptions       const & options,
                      BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    if (options.alphReduction == 10)
    {
        StringSet<String<AminoAcid10>, Owner<ConcatDirect<>>> newQry;
        newQry.concat = oldQry.concat; //implicit conversion
        newQry.limits = oldQry.limits;
        return step04_seedAndSearch(matchManager, newQry,
                                           dbIndex, options, TFormat());
    }

    return step04_seedAndSearch(matchManager, oldQry, dbIndex, options, TFormat());
}

// --------------------------------------------------------------------------
// Function preprocessQuery()
// --------------------------------------------------------------------------

template <typename TNewQry,
          typename TOldQry,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step02a_preprocessQuery(TNewQry & newQry,
                                    TOldQry const & oldQry,
                                    ScarletOptions const & options,
                                    BlastFormat<m,
                                                BlastFormatOptions::BlastN,
                                                g> const & /*tag*/)
{
    if (options.revComp)
        reserve(newQry.concat, 2*length(oldQry.concat), Exact());
    else
        reserve(newQry.concat, length(oldQry.concat), Exact());
    for (auto const & s : oldQry)
    {
        appendValue(newQry, s, Exact());
        if (options.revComp)
        {
            appendValue(newQry, s, Exact());
            reverseComplement(back(newQry));
            //TODO check if the above actually works
        }
    }
}

template <typename TNewQry,
          typename TOldQry,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step02a_preprocessQuery(TNewQry & newQry,
                                    TOldQry const & oldQry,
                                    ScarletOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::BlastP,
                                                g> const & /*tag*/)
{
    newQry.concat = oldQry.concat; //implicit conversion
    newQry.limits = oldQry.limits;
}

template <typename TNewQry,
          typename TOldQry,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step02a_preprocessQuery(TNewQry & newQry,
                                    TOldQry const & oldQry,
                                    ScarletOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::TBlastN,
                                                g> const & /*tag*/)
{
    newQry.concat = oldQry.concat; //implicit conversion
    newQry.limits = oldQry.limits;
}

template <typename TNewQry,
          typename TOldQry,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step02a_preprocessQuery(TNewQry & newQry,
                                    TOldQry const & oldQry,
                                    ScarletOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::BlastX,
                                                g> const & /*tag*/)
{
    reserve(newQry.concat, length(oldQry.concat) / 3 * 6, Exact());
    //  / 3 * 6 is a little too much, but the speedup is worth it
    for (auto const & s : oldQry)
        appendSixFrameTranslation(newQry, s);
}

template <typename TNewQry,
          typename TOldQry,
          BlastFormatOptions::M m,
          BlastFormatOptions::Generation g>
inline void step02a_preprocessQuery(TNewQry & newQry,
                                    TOldQry const & oldQry,
                                    ScarletOptions const & /**/,
                                    BlastFormat<m,
                                                BlastFormatOptions::TBlastX,
                                                g> const & /*tag*/)
{
    reserve(newQry.concat, length(oldQry.concat) / 3 * 6, Exact());
    //  / 3 * 6 is a little too much, but the speedup is worth it
    for (auto const & s : oldQry)
        appendSixFrameTranslation(newQry, s);
}

template <typename TDbIndex,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step02_preprocessQuery(MatchManager               & matchManager,
                       StringSet<CharString, Owner<ConcatDirect<> > >  & oldQry,
                       TDbIndex                   & dbIndex,
                       ScarletOptions       const & options,
                       BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    typename UnreducedStringSet<p>::Type newQry;

    double start = sysTime();
    std::cout << "Preprocessing Query Sequences..." << std::flush;

    // depending on blastProgram
    step02a_preprocessQuery(newQry, oldQry, options, TFormat());

    clear(oldQry);
    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "Number of queries after preproc: " << length(newQry)
                << "\n\n" << std::flush;

    start = sysTime();
    std::cout << "Dumping unreduced Query Sequences..." << std::flush;

    //TODO save to TMPDIR instead
    CharString _path = options.queryFile;
    append(_path, ".unredqry");
    save(newQry, toCString(_path));
    std::cout << " done.\n";
    finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;

    if (options.alphReduction != 0)
        return step03_reduceAlphabet(matchManager, newQry, dbIndex,
                                     options, TFormat());
    else
        return step04_seedAndSearch(matchManager, newQry, dbIndex,
                                    options, TFormat());
}

// --------------------------------------------------------------------------
// Function loadQuerySeuences()
// --------------------------------------------------------------------------


template <typename TDbIndex,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step01_loadQuerySequences(MatchManager & mm,
                          TDbIndex & dbIndex,
                          ScarletOptions const & options,
                          BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    // load sequences as CharString first to be error tolerant
    StringSet<CharString, Owner<ConcatDirect<> > > qrySeqs;
    double start = sysTime();
    std::cout << "Loading Query Sequences..." << std::flush;
    int ret = 0;
    if (options.fileFormat == 1)
        ret = loadSequences(qrySeqs, options.queryFile, Fastq());
    else
        ret = loadSequences(qrySeqs, options.queryFile, Fasta());
    if (ret)
        return ret;
    std::cout << " done.\n";
    double finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;

    unsigned long maxLen = 0ul;
    for (auto const & s : qrySeqs)
        if (length(s) > maxLen)
            maxLen = length(s);
    std::cout << "Number of sequences read: " << length(qrySeqs)
              << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;

    return step02_preprocessQuery(mm, qrySeqs, dbIndex, options, TFormat());
}

// --------------------------------------------------------------------------
// Function loadIndexFromDisk()
// --------------------------------------------------------------------------

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
step00_loadDbIndexFromDisk(MatchManager & mm,
                           ScarletOptions const & options,
                           BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    typedef typename UnreducedStringSet<p>::Type TSubjSeqs;
    typedef Index<TSubjSeqs,      IndexSa<> >      TDbIndex;
    TDbIndex dbIndex;

    std::cout << "Loading Database Index from disk..." << std::flush;
    double start = sysTime();
    int ret = open(dbIndex, toCString(options.dbFile));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }
    double finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    std::cout << "No of Fibres: " << length(indexSA(dbIndex))
              << "\t no of Seqs in Db: " << length(indexText(dbIndex))
              << "\n\n";

    return step01_loadQuerySequences(mm, dbIndex, options, TFormat());
}



// --------------------------------------------------------------------------
// Function searchPipeline()        -- start the search
// --------------------------------------------------------------------------


template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
searchPipeline(MatchManager & mm,
              ScarletOptions const & options,
              BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;
    return step00_loadDbIndexFromDisk(mm, options, TFormat());
}

#endif // header guard
