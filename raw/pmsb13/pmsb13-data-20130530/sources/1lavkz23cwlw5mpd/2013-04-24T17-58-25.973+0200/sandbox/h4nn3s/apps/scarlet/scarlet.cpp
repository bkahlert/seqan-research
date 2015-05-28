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
// scarlet_indexer.cpp: Main File for the indexer application
// ==========================================================================

#undef SEQAN_HAS_ZLIB

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <iostream>

#include "store.hpp"
#include "finder.hpp"


using namespace seqan;


// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class ScarletOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct ScarletOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;


    CharString  queryFile;
    int         queryFormat; // 0 -- fasta, 1 -- fastq
    int         queryAlph; // 0 -- nucleotide, 1 -- aminoacid
    bool        queryTranslate;

    CharString  dbFile;
    int         dbAlph; // 0 -- nucleotide, 1 -- aminoacid
    bool        dbTranslate;

    int         alphReduction;

    CharString  output;
    
    ScarletOptions() :
        verbosity(1)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(ScarletOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("scarlet");
    // Set short description, version, and date.
    setShortDescription(parser, "BLAST compatible local aligner optimized for "
                                "NGS and Metagenomics.");
    setVersion(parser, "0.1");
    setDate(parser, "April 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\f-q QUERY.fasta\\fP "
                         "\\f-d DATABASE.fasta\\fP "
                         "[\\f-o output.m8\\fP]");
    addDescription(parser, "Scarlet is a local aligner that is faster than "
                           "BLAST, optimized for many query sequences.");
    
    addSection(parser, "Input Options");
    addOption(parser, seqan::ArgParseOption("q",
                                            "query",
                                            "Query sequences (fasta).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "q");
    setValidValues(parser, "query", "fasta fa fna faa fastq fq");

    addOption(parser, seqan::ArgParseOption("a",
                                            "alph",
                                            "alphabet type of query sequences.",
                                            seqan::ArgParseArgument::STRING));
    setValidValues(parser, "alph", "auto nucl prot");
    setDefaultValue(parser, "alph", "auto");

    addOption(parser, seqan::ArgParseOption("t",
                                            "translate",
                                            "if query is nucl, translate to "
                                             "protein (six-frame-translation).")
             );
    setDefaultValue(parser, "translate", "false");

    
    addOption(parser, seqan::ArgParseOption("d",
                                            "database",
                                            "Database sequences (fasta), with "
                                             "precomputed index (.sa).",
                                            seqan::ArgParseArgument::INPUTFILE,
                                            "IN"));
    setRequired(parser, "d");
    setValidValues(parser, "database", "fasta fa fna faa");

    // TODO implement some way to transfer options from index builder
    // maybe write a spec file that is loaded by this app

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o",
                                            "output",
                                            "File to hold reports on hits (.m8 "
                                            "is blastall -m8 et cetera)",
                                            seqan::ArgParseArgument::OUTPUTFILE,
                                            "OUT"));
    setValidValues(parser, "output", "m8 m9");
    setDefaultValue(parser, "output", "output.m8");



    addTextSection(parser, "Environment Variables");
    addListItem(parser, "\\fBTMPDIR\\fP",
                        "set this to a local directory with lots of "
                        "space. If you can afford it use /dev/shm.");
    
    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBscarlet\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // Extract option values.
    getOptionValue(options.queryFile, parser, "query");
    getOptionValue(options.dbFile, parser, "database");

    seqan::CharString alph;
    seqan::getOptionValue(alph, parser, "alph");
    //verifyFileFormat(); TODO
    options.dbAlph = 1; // protein

    getOptionValue(options.output, parser, "output");

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
    ScarletOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "SCARLET v0.1\n";

    int ret = 0;

    TQuerySeqs  qrySeqs;
    TDbSeqs     dbSeqs;
    TDbIndex    dbIndex;

    double start = 0;
    double finish = 0;

    // Load Database sequences
    std::cout << "Loading Database sequences from disk..." << std::flush;
    start = sysTime();
    CharString _dbSeqs = options.dbFile;
    append(_dbSeqs, ".txt");
    ret = open(dbSeqs, toCString(_dbSeqs));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }
    finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    unsigned long maxLen = 0ul;
    for (auto const & s : dbSeqs)
        if (length(s) > maxLen)
            maxLen = length(s);
    std::cout << "Number of sequences read: " << length(dbSeqs)
              << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;

    // Load Index
    std::cout << "Loading Database Index from disk..." << std::flush;
    start = sysTime();
    ret = open(dbIndex, toCString(options.dbFile));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }
    finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << finish << "s \n\n" << std::flush;


    // Load Query Sequences
    start = sysTime();
    std::cout << "Loading Query Sequences..." << std::flush;
    if (options.queryFormat == 1)
        ret = loadSequences(qrySeqs, options.queryFile, Fastq());
    else
        ret = loadSequences(qrySeqs, options.queryFile, Fasta());
    if (ret)
        return ret;
    std::cout << " done.\n";
    finish = sysTime() - start;
    std::cout << "Runtime: " << finish << "s \n" << std::flush;
    
    _debug_shorten(qrySeqs, 30);
    maxLen = 0ul;
    for (auto const & s : qrySeqs)
        if (length(s) > maxLen)
            maxLen = length(s);
    std::cout << "Number of sequences read: " << length(qrySeqs)
              << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;


    // Generate Index
    std::cout << "Generating Query-Index..." << std::flush;
    start = sysTime();
    TQueryIndex qryIndex(qrySeqs);
    // we only want full query sequences in index
    typedef typename Fibre<TQueryIndex, EsaSA>::Type TSa;
    TSa & sa = indexSA(qryIndex);
    clear(sa);
    reserve(sa, length(qrySeqs), Exact());
    for (int i = 0; i < length(qrySeqs); ++i)
    {
        Value<TSa>::Type n;
        assignValueI1(n, i);
        assignValueI2(n, length(qrySeqs[i])-1);
        appendValue(sa, n, Exact());
    }
    
    for (auto & sav : indexSA(qryIndex))
        std::cout << "SA pair: " << sav << std::endl << std::flush;
    
    Iterator<TQueryIndex, TopDown<> >::Type it(qryIndex); // instantiate
    finish = sysTime() - start;
    std::cout << " done.\n" << std::flush;
    std::cout << "-------\n";

    for (auto & sav : indexSA(qryIndex))
        std::cout << "SA pair: " << sav << std::endl << std::flush;

    // FIND
    std::cout << "Starting a search...\n" << std::flush;
    MatchManager matchManager;
    ScarletFinder finder(matchManager);

    ext::find(finder, dbIndex, qryIndex, 1);

    for (Match const & m : matchManager.matches)
        _printMatch(m);

    return 0;
}
