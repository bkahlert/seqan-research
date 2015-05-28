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


using namespace seqan;

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

    std::cout << "SCARLET v0.1\n";

    // Load Database sequences
//     std::cout << "Loading Database sequences from disk..." << std::flush;
//     start = sysTime();
//     CharString _dbSeqs = options.dbFile;
//     append(_dbSeqs, ".txt");
//     ret = open(dbSeqs, toCString(_dbSeqs));
//     if (ret != true)
//     {
//         std::cout << " failed.\n" << std::flush;
//         return 1;
//     }
//     finish = sysTime() - start;
//     std::cout << " done.\n" << std::flush;
//     std::cout << "Runtime: " << finish << "s \n" << std::flush;
//     unsigned long maxLen = 0ul;
//     for (auto const & s : dbSeqs)
//         if (length(s) > maxLen)
//             maxLen = length(s);
//     std::cout << "Number of sequences read: " << length(dbSeqs)
//               << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;



    MatchManager matchManager;

    int ret = beginPipeline(matchManager, options);
    if (ret)
        return ret;

    seedMatches2QueryMatches(matchManager, options);

    String<Match> filteredMatches;
    joinAndFilterMatches(filteredMatches, matchManager, options);

    for (Match const & m : filteredMatches)
        _printMatch(m);

    return 0;
}
