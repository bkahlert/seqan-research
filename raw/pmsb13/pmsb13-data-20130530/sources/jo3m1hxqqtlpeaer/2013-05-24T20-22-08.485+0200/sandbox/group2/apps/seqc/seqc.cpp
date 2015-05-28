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
// Author: Daniel Kersting <dkersting@fu-berlin.de>
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
#include "input.h"
#include "output.h"
#include "kmer_content.h"

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------
// Program entry point.
using std::cout;

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

    seqan::SequenceStream seqStreamIn(seqan::toCString(options.inputFile));
    
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
    KmerContent kmerContent(1,3);

    unsigned noReads = 0;

    while((noReads < options.maxReads) && !atEnd(seqStreamIn))
    {
        ++noReads;
        if (readRecord(id, seq, qual, seqStreamIn) == 0){
           readStats.collectReadStats(id, seq, qual); 
	   kmerContent.CountKmers(seq);
           readStats.dupCheck(seq);
        } else {
            std::cerr << "WARN: Bad Sequence in record [#" << noReads << "] " << id << std::endl;
        }
    }



    
    readStats.jobParams["Filename"] = seqan::toCString(options.inputFile); 
    readStats.jobParams["report output dir"] = seqan::toCString(options.outdir); 

    TSVWriter writer;

   
    writer.writeAll(readStats, kmerContent, options.outdir);
    
    std::cout << "no of reads analyzed: " << noReads << std::endl;
    return 0;
}
