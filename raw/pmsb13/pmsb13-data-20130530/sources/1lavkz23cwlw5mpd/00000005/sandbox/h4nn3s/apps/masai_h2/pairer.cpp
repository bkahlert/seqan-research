// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the masai_h2_output_pe application.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "options.h"
#include "pairer.h"
#include "writer.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options : public MasaiOptions
{
    CharString  genomeFile;

    CharString  readsLeftFile;
    CharString  readsRightFile;
    CharString  mappedReadsLeftFile;
    CharString  mappedReadsRightFile;

    CharString  mappedPairsFile;
    OutputFormat outputFormat;
    bool        outputCigar;

//    unsigned    errorsPerRead;
    bool        mismatchesOnly;
    unsigned    libraryLength;
    unsigned    libraryError;

    bool        noDump;

    Options() :
        MasaiOptions(),
        outputFormat(SAM),
        outputCigar(true),
//        errorsPerRead(5),
        mismatchesOnly(false),
        libraryLength(220),
        libraryError(50),
        noDump(true)
    {}
};

// ----------------------------------------------------------------------------
// Function setupArgumentParser()                              [ArgumentParser]
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "masai_h2_output_pe");
    setShortDescription(parser, "Masai Output - Paired End Mode");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE L\\fP> <\\fIREADS FILE R\\fP> <\\fIRAW FILE L\\fP> <\\fIRAW FILE R\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");
    setValidValues(parser, 1, "fastq fasta fa");
    setValidValues(parser, 2, "fastq fasta fa");
    setValidValues(parser, 3, "raw");
    setValidValues(parser, 4, "raw");


    addSection(parser, "Pairing Options");

//    addOption(parser, ArgParseOption("e",  "errors", "Maximum number of errors per read.", ArgParseOption::INTEGER));
//    setMinValue(parser, "errors", "0");
//    setMaxValue(parser, "errors", "32");
//    setDefaultValue(parser, "errors", options.errorsPerRead);

    addOption(parser, ArgParseOption("ng", "no-gaps", "Do not align reads with gaps."));

    addOption(parser, ArgParseOption("ll", "library-length", "Library length.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "library-length", options.libraryLength);

    addOption(parser, ArgParseOption("le", "library-error", "Library length tolerance.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "library-error", options.libraryError);


    addSection(parser, "Output Options");

    setTmpFolder(parser);
    setOutputFile(parser, options);
    addOption(parser, ArgParseOption("nc", "no-cigar", "Do not output CIGAR string. This only affects SAM output."));


    addSection(parser, "Debug Options");

    addOption(parser, ArgParseOption("nd", "no-dump", "Do not dump results."));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()                                        [Options]
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Parse genome input file.
    getArgumentValue(options.genomeFile, parser, 0);

    // Parse reads input files.
    getArgumentValue(options.readsLeftFile, parser, 1);
    getArgumentValue(options.readsRightFile, parser, 2);

    // Parse raw input files.
    getArgumentValue(options.mappedReadsLeftFile, parser, 3);
    getArgumentValue(options.mappedReadsRightFile, parser, 4);

    // Parse pairing options.
//    getOptionValue(options.errorsPerRead, parser, "errors");
    options.mismatchesOnly = isSet(parser, "no-gaps");
    getOptionValue(options.libraryLength, parser, "library-length");
    getOptionValue(options.libraryError, parser, "library-error");

    // Parse tmp folder.
    getTmpFolder(options, parser);

    // Parse output file.
    getOutputFile(options.mappedPairsFile, options, parser, options.readsLeftFile, "_pe");

    // Parse output format.
    getOutputFormat(options, options.mappedPairsFile);
    options.outputCigar = !isSet(parser, "no-cigar");

    // Parse debug options.
    options.noDump = isSet(parser, "no-dump");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function runPairer()
// ----------------------------------------------------------------------------

template <typename TDistance, typename TFormat, typename TReadsConfig>
int runPairer(Options & options, TReadsConfig const & /* config */)
{
    typedef Genome<void>                                                        TGenome;
    typedef Reads<PairedEnd, TReadsConfig>                                      TReads;
    typedef ReadsLoader<PairedEnd, TReadsConfig>                                TReadsLoader;
    typedef Writer<TGenome, TReads, TFormat, TDistance>                         TWriter;
    typedef Pairer<TReads, TWriter>                                             TPairer;

    TFragmentStore      store;
    TGenome             genome(store);
    TReads              reads(store);
    TReadsLoader        readsLoader(reads);
    TWriter             writer(genome, options.noDump);
    TPairer             pairer(writer, options.libraryLength, options.libraryError);

    double start, finish;

    // Load genome.
    std::cout << "Loading genome:\t\t\t" << std::flush;
    start = sysTime();
    if (!load(genome, options.genomeFile))
    {
        std::cerr << "Error while loading genome" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    // Open reads file.
    start = sysTime();
    if (!open(readsLoader, options.readsLeftFile, options.readsRightFile))
    {
        std::cerr << "Error while opening reads file" << std::endl;
        return 1;
    }

    // Open matches files.
    if (!open(pairer, options.mappedReadsLeftFile, options.mappedReadsRightFile))
    {
        std::cerr << "Error while opening matches file" << std::endl;
        return 1;
    }

    // Open output file.
    if (!open(writer, options.mappedPairsFile))
    {
        std::cerr << "Error while opening output file" << std::endl;
        return 1;
    }

    // Configure writer.
    writeAlignments(writer, options.outputCigar);

    // Reserve space for reads.
    reserve(reads);

    // Load all reads.
    std::cout << "Loading reads:\t\t\t" << std::flush;
    start = sysTime();
    load(readsLoader);
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;
    std::cout << "Reads count:\t\t\t" << reads.readsCount << std::endl;

    // Pass reads to writer.
    setReads(writer, reads);

    // Pass reads to pairer.
    setReads(pairer, reads);

    // Pair matches.
    std::cout << "Pairing matches:\t\t" << std::flush;
    start = sysTime();
    pair(pairer);
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;
    std::cout << "Pairs count:\t\t\t" << pairer.pairsCount << std::endl;

    // Close output file.
    close(writer);

    // Close matches file.
    close(pairer);

    // Close reads file.
    close(readsLoader);

    return 0;
}

template <typename TDistance, typename TFormat>
int runPairer(Options & options)
{
    typedef typename IsSameType<TFormat, Sam>::Type                         TUseReadStore;
    typedef typename IsSameType<TFormat, Sam>::Type                         TUseReadNameStore;
    typedef ReadsConfig<TUseReadStore, TUseReadNameStore>                   TReadsConfig;

    return runPairer<TDistance, TFormat>(options, TReadsConfig());
}

// ----------------------------------------------------------------------------
// Functions configure*()
// ----------------------------------------------------------------------------

template <typename TDistance>
int configureOutputFormat(Options & options)
{
    switch (options.outputFormat)
    {
    case Options::RAW:
        return runPairer<TDistance, Raw>(options);

    case Options::SAM:
        return runPairer<TDistance, Sam>(options);

    default:
        return 1;
    }
}

int configureDistance(Options & options)
{
    if (options.mismatchesOnly)
        return configureOutputFormat<HammingDistance>(options);
    else
        return configureOutputFormat<EditDistance>(options);
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    return configureDistance(options);
}
