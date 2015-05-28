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
// scarlet.cpp: Main File for Scarlet
// ==========================================================================

#undef SEQAN_HAS_ZLIB

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <iostream>

#include "options.hpp"
#include "scarlet.hpp"
#include "finder.hpp"
#include "writer.hpp"


using namespace seqan;


// --------------------------------------------------------------------------
// Function mainComp()
// --------------------------------------------------------------------------

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
mainComp(ScarletOptions const & options,
         BlastFormat<m,p,g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;

    std::cout << "SCARLET v0.1\n";

    String<Match> filteredMatches;

    {
        MatchManager matchManager;

        // find matches
        int ret = searchPipeline(matchManager, options, TFormat());
        if (ret)
            return ret;

        // post-process matches
        seedMatches2QueryMatches(matchManager, options);
        joinAndFilterMatches(filteredMatches, matchManager, options);

    } // destruct matchManager with members


    for (Match const & ma : filteredMatches)
        _printMatch(ma);

//     return 0;
    return writePipeline(filteredMatches, options, TFormat());
}


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    ScarletOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // CONVERT Run-time options to compile-time Format-Type
    switch (options.blastProg)
    {
        case BlastFormatOptions::BlastN :
            if (hasSuffix(options.output, ".m8"))
            {
                typedef BlastFormat<BlastFormatOptions::Tabular,
                                    BlastFormatOptions::BlastN,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            } else if (hasSuffix(options.output, ".m9"))
            {
                typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                                    BlastFormatOptions::BlastN,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            }
            break;
        case BlastFormatOptions::BlastP :
            if (hasSuffix(options.output, ".m8"))
            {
                typedef BlastFormat<BlastFormatOptions::Tabular,
                                    BlastFormatOptions::BlastP,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            } else if (hasSuffix(options.output, ".m9"))
            {
                typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                                    BlastFormatOptions::BlastP,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            }
            break;
        case BlastFormatOptions::BlastX :
            if (hasSuffix(options.output, ".m8"))
            {
                typedef BlastFormat<BlastFormatOptions::Tabular,
                                    BlastFormatOptions::BlastX,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            } else if (hasSuffix(options.output, ".m9"))
            {
                typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                                    BlastFormatOptions::BlastX,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            }
            break;
        case BlastFormatOptions::TBlastN :
            if (hasSuffix(options.output, ".m8"))
            {
                typedef BlastFormat<BlastFormatOptions::Tabular,
                                    BlastFormatOptions::TBlastN,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            } else if (hasSuffix(options.output, ".m9"))
            {
                typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                                    BlastFormatOptions::TBlastN,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            }
            break;
        case BlastFormatOptions::TBlastX :
            if (hasSuffix(options.output, ".m8"))
            {
                typedef BlastFormat<BlastFormatOptions::Tabular,
                                    BlastFormatOptions::TBlastX,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            } else if (hasSuffix(options.output, ".m9"))
            {
                typedef BlastFormat<BlastFormatOptions::TabularWithHeader,
                                    BlastFormatOptions::TBlastX,
                                    BlastFormatOptions::Blast> format;
                return mainComp(options, format());
            }
            break;
        default:
            return -1;
    }
    return -1;
}
