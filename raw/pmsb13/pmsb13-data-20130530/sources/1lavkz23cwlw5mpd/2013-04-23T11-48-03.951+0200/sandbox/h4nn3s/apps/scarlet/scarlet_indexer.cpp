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
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/index_extras.h>

#include <iostream>

#include "store.hpp"


using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct ScarletIndexerOptions
{
//     // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
//     int verbosity;

    seqan::CharString   dbFile;
    int                 dbAlph; // 0 -- nucleotide, 1 -- aminoacid

    seqan::CharString   indexFile;

    bool                translate;


    ScarletIndexerOptions()
    {
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(ScarletIndexerOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("scarlet_indexer");
    // Set short description, version, and date.
    setShortDescription(parser, "TODO Put a Short Description Here");
    setVersion(parser, "TODO0.1");
    setDate(parser, "TODOJuly 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\-i DATABASE.fasta\\fP [\\-o INDEXFILE.sa\\fP]");
    addDescription(parser, "This is the application skelleton and you should modify this string.");


    addSection(parser, "Input Options");
    addOption(parser, seqan::ArgParseOption("i",
                                            "input",
                                            "Database sequences (fasta).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "i");
    setValidValues(parser, "input", "fasta fa fna faa");

    addOption(parser, seqan::ArgParseOption("a",
                                            "alph",
                                            "alphabet type of database sequences",
                                            seqan::ArgParseArgument::STRING));
    setValidValues(parser, "alph", "auto nucl prot");
    setDefaultValue(parser, "alph", "auto");

    addOption(parser, seqan::ArgParseOption("t",
                                            "translate",
                                            "if input is nucl, translate to "
                                             "protein (six-frame-translation).")
             );
    setDefaultValue(parser, "translate", "false");


    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o",
                                            "output",
                                            "Index of database sequences",
                                            seqan::ArgParseArgument::OUTPUTFILE,
                                            "OUT"));
    setValidValues(parser, "output", "sa fm");



    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");

    addTextSection(parser, "Remarks");
    addText(parser, "Note that the indeces created are binary and not "
                    "necessarily platform independant.");
    addText(parser, "Note that the indeces created are binary and not "
                    "necessarily platform independant.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);


    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // Extract option values.

    seqan::getOptionValue(options.dbFile, parser, "input");

    seqan::CharString alph;
    seqan::getOptionValue(alph, parser, "alph");
    //verifyFileFormat(); TODO
    options.dbAlph = 1; // protein

    if (seqan::isSet(parser, "output"))
        seqan::getOptionValue(options.indexFile, parser, "output");
    else
    {
        options.indexFile = options.dbFile;// appends ending automatically
    }

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

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

    std::cout << "Scarlet Indexer\n"
              << "===============\n\n";


    double s = 0;
    double e = 0;

    // Load Sequences
    s = sysTime();
    std::cout << "Loading Sequences..." << std::flush;
    TDbSeqs seqs;
    int re = loadSequences(seqs, options.dbFile, Fasta());
    if (re)
        return re;
    std::cout << " done.\n";
    e = sysTime() - s;
    std::cout << "Runtime: " << e << "s \n" << std::flush;
    unsigned long maxLen = 0ul;
    for (auto it = begin(seqs), itEnd = end(seqs); it != itEnd; ++it)
        if (length(*it) > maxLen)
            maxLen = length(*it);
    std::cout << "Number of sequences read: " << length(seqs)
              << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;


    // Generate Index
    std::cout << "Generating Index..." << std::flush;
    s = sysTime();
    TDbIndex dbIndex(seqs);
    Iterator<TDbIndex, TopDown<> >::Type it(dbIndex); // instantiate
    e = sysTime() - s;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << e << "s \n\n" << std::flush;


    // Dump Index
    std::cout << "Writing Index to disk..." << std::flush;
    s = sysTime();
    save(dbIndex, toCString(options.indexFile));
    e = sysTime() - s;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << e << "s \n" << std::flush;

    return 0;
}
