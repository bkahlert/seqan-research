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


struct SharedOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;


    CharString  dbFile;

    // for indexer, the file format of database sequences
    // for main app, the file format of query sequences
    // 0 -- fasta, 1 -- fastq
    int         fileFormat; 

    int         alphReduction;

    BlastFormatOptions::Program blastProg;

    SharedOptions() :
        verbosity(1), alphReduction(0), blastProg(BlastFormatOptions::BlastX),
        fileFormat(0)
    {}
};


struct ScarletOptions : public SharedOptions
{

    CharString  queryFile;
    bool        revComp;

    CharString  output;


    unsigned char seedLength;
    unsigned char maxSeedDist;

    bool        semiGlobal;


    ScarletOptions() :
        SharedOptions(), seedLength(30), maxSeedDist(1), revComp(false),
        semiGlobal(false)
    {}
};

struct ScarletIndexerOptions : public SharedOptions
{
    ScarletIndexerOptions()
        : SharedOptions()
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


    //TODO
    options.blastProg = BlastFormatOptions::BlastP;

    getOptionValue(options.output, parser, "output");

    return seqan::ArgumentParser::PARSE_OK;
}


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


    //verifyFileFormat(); TODO



    return seqan::ArgumentParser::PARSE_OK;
}


#endif // header guard