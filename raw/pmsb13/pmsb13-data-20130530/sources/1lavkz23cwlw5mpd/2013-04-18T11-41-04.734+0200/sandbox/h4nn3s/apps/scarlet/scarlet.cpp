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
#include <seqan/index.h>
#include <seqan/index_extras.h>

#include <iostream>


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

template <typename TString>
inline int
loadSequences(StringSet<TString, Owner<ConcatDirect<> > > & seqs,
              ScarletOptions const & options)
{
    std::ifstream stream;
    stream.open(toCString(options.queryFile));

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    StringSet<CharString, Owner<ConcatDirect<> > > ids;

    int res = read2(ids, seqs, reader, Fasta());
    if (res)
    {
        std::cerr << "Error : " << res << "\n";
        return res;
    }

    stream.close();
    return 0;
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

    typedef StringSet<CharString,  Owner<ConcatDirect<> > > TSeqs;
    typedef Index<TSeqs, IndexSa<> > TIndex;

    TSeqs   qrySeqs;
    TSeqs   dbSeqs;
    TIndex  dbIndex;

    double s = 0;
    double e = 0;

    // Load Database sequences
    std::cout << "Loading Database sequences from disk..." << std::flush;
    s = sysTime();
    CharString _dbSeqs = options.dbFile;
    append(_dbSeqs, ".txt");
    ret = open(dbSeqs, toCString(_dbSeqs));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }
    e = sysTime() - s;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << e << "s \n" << std::flush;
    unsigned long maxLen = 0ul;
    for (auto it = begin(dbSeqs), itEnd = end(dbSeqs); it != itEnd; ++it)
        if (length(*it) > maxLen)
            maxLen = length(*it);
    std::cout << "Number of sequences read: " << length(dbSeqs)
              << "\nLongest sequence read: " << maxLen << "\n\n" << std::flush;

    // Load Index
    std::cout << "Loading Database Index from disk..." << std::flush;
    s = sysTime();
    ret = open(dbIndex, toCString(options.dbFile));
    if (ret != true)
    {
        std::cout << " failed.\n" << std::flush;
        return 1;
    }
    e = sysTime() - s;
    std::cout << " done.\n" << std::flush;
    std::cout << "Runtime: " << e << "s \n" << std::flush;


    return 0;
}
