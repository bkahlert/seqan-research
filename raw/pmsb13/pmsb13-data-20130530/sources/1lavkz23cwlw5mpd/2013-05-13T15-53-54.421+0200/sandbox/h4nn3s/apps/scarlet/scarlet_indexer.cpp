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
// scarlet.cpp: Main File for the main application
// ==========================================================================

// why is this neccessary?
#undef SEQAN_HAS_ZLIB

// #define SEQAN_DEBUG_INDEX

#include <seqan/basic.h>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>

#include "scarlet_indexer.hpp"

using namespace seqan;

// ==========================================================================
// Functions
// ==========================================================================
// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
int mainComp(ScarletIndexerOptions const & options,
             BlastFormat<m,p,g> const & /*tag*/)
{
    typedef BlastFormat<m,p,g> TFormat;
    std::cout << "Scarlet Indexer\n"
              << "===============\n\n";

    return beginPipeline(options, TFormat());
}


// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    ScarletIndexerOptions options;
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
            {
                typedef BlastFormat<BlastFormatOptions::INVALID_M,
                                    BlastFormatOptions::BlastN,
                                    BlastFormatOptions::INVALID_Generation> format;
                return mainComp(options, format());
            }break;
        case BlastFormatOptions::BlastP :
            {
                typedef BlastFormat<BlastFormatOptions::INVALID_M,
                                    BlastFormatOptions::BlastP,
                                    BlastFormatOptions::INVALID_Generation> format;
                return mainComp(options, format());
            }
            break;
        case BlastFormatOptions::BlastX :
            {
                typedef BlastFormat<BlastFormatOptions::INVALID_M,
                                    BlastFormatOptions::BlastX,
                                    BlastFormatOptions::INVALID_Generation> format;
                return mainComp(options, format());
            }
            break;
        case BlastFormatOptions::TBlastN :
            {
                typedef BlastFormat<BlastFormatOptions::INVALID_M,
                                    BlastFormatOptions::TBlastN,
                                    BlastFormatOptions::INVALID_Generation> format;
                return mainComp(options, format());
            }
            break;
        case BlastFormatOptions::TBlastX :
            {
                typedef BlastFormat<BlastFormatOptions::INVALID_M,
                                    BlastFormatOptions::TBlastX,
                                    BlastFormatOptions::INVALID_Generation> format;
                return mainComp(options, format());
            }
            break;
        default:
            return -1;
    }
    return -1;
}
