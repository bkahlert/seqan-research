// ==========================================================================
//                                   SeqDPT
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Sebastian Roskosch, Benjamin Strauch
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------
// Program entry point.

int main(int argc, char const ** argv)
{
	// PARSER DEFINITON ----------------------------------------
	seqan::ArgumentParser parser("SeqDPT");
	setDate(parser, __DATE__);
	addArgument(parser, seqan::ArgParseArgument(
			seqan::ArgParseArgument::INPUTFILE, "FILES", true));

	// READ TRIMMING
	addSection(parser, "Trimming options");
	seqan::ArgParseOption qualOpt = seqan::ArgParseOption(
			"q", "quality", "Quality threshold for read trimming",
			seqan::ArgParseArgument::INTEGER, "PHRED");
	setRequired(qualOpt, true);
	setMinValue(qualOpt, "0");
	setMaxValue(qualOpt, "40");
	addOption(parser, qualOpt);

	seqan::ArgParseOption lenOpt = seqan::ArgParseOption(
			"l", "length", "Minimum read length after trimming. Shorter reads will be removed.",
			seqan::ArgParseArgument::INTEGER, "LENGTH");
	setDefaultValue(lenOpt, 1);
	addOption(parser, lenOpt);

	seqan::ArgParseOption trimOpt = seqan::ArgParseOption(
			"m", "method", "Method for trimming reads.",
			seqan::ArgParseArgument::STRING, "METHOD");
	setDefaultValue(trimOpt, "WIN");
	setValidValues(trimOpt, "WIN BWA TAIL");
	addOption(parser, trimOpt);
	// END PARSER DEFINITON ----------------------------------------

	// Additional checks
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// Check if input was successfully parsed.
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	// Check if one or two input files (single or paired-end) were given.
	int fileCount = getArgumentValueCount(parser, 0);
	if (!(fileCount == 1 || fileCount == 2)){
		printShortHelp(parser);
		return seqan::ArgumentParser::PARSE_HELP;
	}

	// TODO: Different procedures for handling one file or two files.
	seqan::CharString fileName1;
	int cutoff, min_length;

	getArgumentValue(fileName1, parser, 0, 0);
	getOptionValue(cutoff, parser, "q");
	getOptionValue(min_length, parser, "l");

    return 0;
}
