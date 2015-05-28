// ==========================================================================
//                                    seqc
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Daniel Kersting <@fu-berlin.de>
// Author: Antje Oldenburg <antje.oldenburg@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
 

#include <seqan/seq_io.h>

#include "seqc.h"
#include "read_stats.h"
#include "output.h"

using namespace seqan;
using std::cout;
using std::endl;
using std::ostream;

void AppOptions::out(ostream & os){
    os << "__ OPTIONS ____________________________________________________________________" << endl
       << "VERBOSITY\t" << verbosity << endl
       << "TEXT     \t" << text << endl;
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("seqc");
    // Set short description, version, and date.
    setShortDescription(parser, "Quality Control app for NGS reads - implemented using seqan lib.");
    setVersion(parser, "0.1");
    setDate(parser, "May 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application description.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
    
    // custom options
    addOption(parser, seqan::ArgParseOption(
                "r","max-reads",
                "Maximum number of records to analyse. Report on what we found so far, as if the input ended after given number of reads",  
                seqan::ArgParseArgument::INTEGER, "INT")); 

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBseqc\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
    seqan::getArgumentValue(options.text, parser, 0);

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
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
            options.out(cout); 
    

    // start reading the input stream
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString qual;

    seqan::SequenceStream seqStreamIn(argv[1]);
    
    if (!isGood(seqStreamIn))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
   

    /*
     * Initialize the data structure for stats 
     * 
     */
  
    ReadStats readStats(1);

    unsigned noReads = 0;
    int writeResult = 0;
    while((noReads < options.maxReads) && !atEnd(seqStreamIn))
    {
        ++noReads;
        if (readRecord(id, seq, qual, seqStreamIn) == 0){
           readStats.collectReadStats(id, seq, qual); 
           readStats.dupCheck(seq);
        } else {
            std::cerr << "ERROR: Bad Sequence" << std::endl;
        }
    }

    
    readStats.jobParams["Filename"] = argv[1]; 
    readStats.jobParams["File type"] = "Conventional base calls";
    readStats.jobParams["Encoding"] = "Sanger / Illumina 1.9";

    TSVWriter writer;

    const std::string outdirname = "testout";
    writer.writeAll(readStats, outdirname);
    
    std::cout << "no of reads analyzed: " << noReads << std::endl;
    return 0;
}
