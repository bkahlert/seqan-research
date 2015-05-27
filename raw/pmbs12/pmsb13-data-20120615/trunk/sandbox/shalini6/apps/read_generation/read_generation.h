// ==========================================================================
//                               readGeneration
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef SANDBOX_SINGER_APPS_READGENERATION_READGENERATION_H_
#define SANDBOX_SINGER_APPS_READGENERATION_READGENERATION_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/basic.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/random.h>
#include <seqan/stream.h>

const int SEED = 42;

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Options
{
    bool showHelp;
    bool showVersion;
    CharString inputFileName;
    CharString outputFileName;
    unsigned readLength;
    unsigned coverage;
    unsigned numRepeats;
    unsigned repeatLength;
    String<CharString> texts;
    
    Options() :
        showHelp(false),
        showVersion(false),
        inputFileName(),
        outputFileName(),
        readLength(32),
        coverage(30),
        numRepeats(5),
        repeatLength(50),
        texts()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

void
setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
    addVersionLine(parser, "0.1");
    
    addTitleLine(parser, "**********************");
    addTitleLine(parser, "* readGeneration *");
    addTitleLine(parser, "**********************");
    addTitleLine(parser, "");
    addTitleLine(parser, "(c) 2012 by Your Name <your.email@example.net>");

    addUsageLine(parser, "[OPTIONS] TEXT+");
    
	addSection(parser, "Main Options");
	addOption(parser, CommandLineOption("if",  "input-file",  "Input file.", OptionType::String | OptionType::Label, options.inputFileName));
	addOption(parser, CommandLineOption("of",  "output-file", "Output file.", OptionType::String | OptionType::Label, options.outputFileName));
	addOption(parser, CommandLineOption("rl",  "read-length", "Read length.", OptionType::Integer | OptionType::Label, options.readLength));
	addOption(parser, CommandLineOption("c",  "coverage", "Coverage.", OptionType::Integer | OptionType::Label, options.coverage));
	addOption(parser, CommandLineOption("n",  "num-repeats", "Num repeats.", OptionType::Integer | OptionType::Label, options.numRepeats));
	addOption(parser, CommandLineOption("r",  "repeat-length", "Repeat length.", OptionType::Integer | OptionType::Label, options.repeatLength));
    //requiredArguments(parser, 1);
}

int parseCommandLineAndCheck(Options & options,
                             CommandLineParser & parser,
                             int argc,
                             char const ** argv)
{
    bool stop = !parse(parser, argc, argv);

    std::cerr << "TEST " << stop << std::endl;

    if (stop)
        return 1;
    if (isSetLong(parser, "help"))
    {
        options.showHelp = true;
        return 0;
    }
    if (isSetLong(parser, "version"))
    {
        options.showVersion = true;
        return 0;
    }

    getOptionValueLong(parser, "input-file", options.inputFileName);
    getOptionValueLong(parser, "output-file", options.outputFileName);
    getOptionValueLong(parser, "read-length", options.readLength);
    getOptionValueLong(parser, "coverage", options.coverage);
    getOptionValueLong(parser, "num-repeats", options.numRepeats);
    getOptionValueLong(parser, "repeat-length", options.repeatLength);
    
    options.texts = getArgumentValues(parser);

	return 0;
}

int mainWithOptions(Options & options)
{
    typedef Dna5 TChar;
    typedef String<TChar> TString;
    typedef Iterator<String<CharString> >::Type TIterator;

    std::cout << "Option Arguments:" << std::endl;
    std::cout << "  input file:  \"" << options.inputFileName << "\"" << std::endl;
    std::cout << "  output file: \"" << options.outputFileName << "\"" << std::endl;
    std::cout << "Non-option Arguments:" << std::endl;

    for (TIterator it = begin(options.texts); it != end(options.texts); ++it)
    {
        std::cout << "  " << *it << std::endl;
    }

    ::std::fstream fstrm;
	fstrm.open(toCString(options.inputFileName), ::std::ios_base::in | ::std::ios_base::binary);
    TString genome;

    read(fstrm, genome, Fasta());
	fstrm.close();

    std::cerr << "genome length: " << length(genome) << std::endl;

	Rng<MersenneTwister> rng(SEED);
	unsigned genomeLength = length(genome);
    StringSet<TString> reads;
    StringSet<String<char> > readIds;
	
	// genome creation
    for (unsigned i = 0; i < options.numRepeats; ++i)
    {
        unsigned repeatPos = pickRandomNumber(rng) % (genomeLength - options.repeatLength - 1);
        TString repeat = infix(genome, repeatPos, repeatPos +  options.repeatLength);
        //std::cerr << "repeat length: " << length(repeat) << " repeatPos: " << repeatPos << " repeatEnd: " <<  repeatPos +  options.repeatLength << std::endl; 
        for (unsigned j = 0; j <= i; ++j)
        {
            repeatPos = pickRandomNumber(rng) % (genomeLength - options.repeatLength - 1);
            //std::cerr << "repeat length: " << length(repeat) << " repeatPos: " << repeatPos << " repeatEnd: " <<  repeatPos +  options.repeatLength << std::endl; 
            insert(genome, repeatPos, repeat, Generous());
        }
    }
    
    std::cerr << "genome length: " << length(genome) << std::endl;

	// read creation
    unsigned readNum = options.coverage * length(genome) / options.readLength;
    for (unsigned i = 0; i < readNum; ++i)
    {
        unsigned readPos = pickRandomNumber(rng) % (genomeLength - options.readLength - 1);
        TString read = infix(genome, readPos, readPos + options.readLength);
        appendValue(reads, read);
        std::stringstream str;
        str << i;
        appendValue(readIds, str.str());
    }

    {
        ::std::FILE * fl = ::std::fopen(toCString(options.outputFileName), "wb");
        write(fl, genome, "a test file", Fasta());
        close (fl);
    }
    {
        ::std::ofstream  fl("readstest.fasta", std::ios_base::out);
        Stream<CharArray< char *> > stream;
        //open(stream, "reads.fasta", OPEN_RDWR | OPEN_CREATE | OPEN_APPEND);

        write2(fl, readIds, reads, Fasta());
        fl.close();
    }
    return 0;
}

#endif  // #ifndef SANDBOX_SINGER_APPS_READGENERATION_READGENERATION_H_

