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
// options.h: contains the options and argument parser
// ==========================================================================


#ifndef SEQAN_SCARLET_OPTIONS_H_
#define SEQAN_SCARLET_OPTIONS_H_

#include <seqan/basic.h>
#include <seqan/blast.h>
#include <seqan/arg_parse.h>

#include <seqan/blast.h>

using namespace seqan;



// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class ScarletOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.

struct ScarletOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;


    CharString  queryFile;
    int         queryFormat; // 0 -- fasta, 1 -- fastq

    bool        revComp; // create reverse complements (blastn-only)

    CharString  dbFile;
    int         dbAlph; // 0 -- nucleotide, 1 -- aminoacid

    int         alphReduction;

    CharString  output;


    unsigned char seedLength;
    unsigned char maxSeedDist;

    BlastFormatOptions::Program blastProg;

    ScarletOptions() :
        verbosity(1), seedLength(30), maxSeedDist(1)
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

    // TODO
    options.queryFormat = 0; // fasta

    getOptionValue(options.output, parser, "output");

    return seqan::ArgumentParser::PARSE_OK;
}




#endif // header guard