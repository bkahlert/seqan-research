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
// Author: Antje Oldenburg <antje.oldenburg@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

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

struct AppOptions
{
        // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The first (and only) argument of the program is stored here.
    CharString text;

    AppOptions() :
        verbosity(1)
    {}
};

struct ReadStats
{

    // 2-dim [int:pos][int:score] array to store counts for
    // scores for each position to calculate 
    // distribution metrics later
    String<String <unsigned> > scoreCount;
    


    ReadStats(unsigned readLength)
    {
        resize(scoreCount, readLength);
        resize(scoreCount[0], 255, 0); 
        std::cout << "Init readStats with readLength " << readLength << std::endl;
    }

    void toString(){
        std::cout << "== Read Stats == " << std::endl;
        for(unsigned i = 0; i < length(scoreCount); ++i)
            std::cout << scoreCount[i] << std::endl;
    }
};


// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("seqc");
    // Set short description, version, and date.
    setShortDescription(parser, "Seqan Quality Control app for NGS reads.");
    setVersion(parser, "0.1");
    setDate(parser, "May 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

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

    std::cout << "Input: " << options.text
                
              << "=============== \n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
            std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "TEXT     \t" << options.text << "\n\n";
    }


    /*
       @TODO
        take input names from args to queue them
    */
    

    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString qual;

    seqan::SequenceStream seqStreamIn(argv[1]);
    seqan::SequenceStream seqStreamOut("compressed.fq.gz", seqan::SequenceStream::WRITE);
    
    if (!isGood(seqStreamIn))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
   
    /*
     * We can guess the scoring system ... @TODO find out which are which
     * Illumina: offset=64  [ 0 - 40] -> [64-104]
     * Sanger:   offset=33  [ 0 - 93] -> [33-126]
     * Solexa:   offset=64  [-5 - 62] -> [59-99]
     *
     */
    unsigned qscoreOffset = 0;

    /*
     * Initialize the data structure for stats 
     * 
     */
  
    ReadStats readStats(999);

    unsigned maxNoReads = 10;
    unsigned noReads = 0;
    int writeResult = 0;
    while((noReads < maxNoReads) && !atEnd(seqStreamIn))
    {
        ++noReads;
        if (readRecord(id, seq, qual, seqStreamIn) == 0){
            std::cout << id << '\t' << seq << '\t' << qual << '\n';
            /*
             * This is where we do the work on each read
             *
             */
            

            /*
             *
             * This is where we decide if we want to filter our output
             */
            writeResult = writeRecord(seqStreamOut, id, seq, qual);
            

        } else {
    /*
     * @TODO what should we do with invalid base calls or/quality values?
     */
            std::cerr << "ERROR: Bad Sequence" << std::endl;
        }
    }
    readStats.toString();
    std::cout << "no of reads analyzed: " << noReads << std::endl;
    return 0;
}
