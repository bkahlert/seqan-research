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
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

#include "readTrimming.h"
#include "adapterTrimming.h"

#include "demultiplex.h"

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------
// Program entry point.

int main(int argc, char const ** argv)
{
	// PARSER DEFINITON ----------------------------------------

	// GENERAL OPTIONS
	seqan::ArgumentParser parser("SeqDPT");
	setDate(parser, __DATE__);
	addArgument(parser, seqan::ArgParseArgument(
			seqan::ArgParseArgument::INPUTFILE, "FILES", true));

	// Barcode Demultiplexing
	addSection(parser, "Demultiplexing Options");
	addOption(parser, seqan::ArgParseOption(
        "app", "approximate", "Select approximate barcodes demultiplexing, allowing one mismatch."));

	seqan::ArgParseOption lenOpt = seqan::ArgParseOption(
			"l", "length", "Minimum read length after trimming. Shorter reads will be removed.",
			seqan::ArgParseArgument::INTEGER, "LENGTH");
	setDefaultValue(lenOpt, 1);
	setMinValue(lenOpt, "1");
	addOption(parser, lenOpt);

	seqan::ArgParseOption trimOpt = seqan::ArgParseOption(
			"m", "method", "Method for trimming reads.",
			seqan::ArgParseArgument::STRING, "METHOD");
	setDefaultValue(trimOpt, "WIN");
	setValidValues(trimOpt, "WIN BWA TAIL");
	addOption(parser, trimOpt);

	// READ TRIMMING
	addSection(parser, "Quality trimming options");
	seqan::ArgParseOption qualOpt = seqan::ArgParseOption(
			"q", "quality", "Quality threshold for read trimming.",
			seqan::ArgParseArgument::INTEGER, "PHRED");
	setRequired(qualOpt, true);
	setMinValue(qualOpt, "0");
	setMaxValue(qualOpt, "40");
	addOption(parser, qualOpt);

	seqan::ArgParseOption lenOpt = seqan::ArgParseOption(
			"l", "length", "Minimum read length after trimming. Shorter reads will be removed.",
			seqan::ArgParseArgument::INTEGER, "LENGTH");
	setDefaultValue(lenOpt, 1);
	setMinValue(lenOpt, "1");
	addOption(parser, lenOpt);

	seqan::ArgParseOption trimOpt = seqan::ArgParseOption(
			"m", "method", "Method for trimming reads.",
			seqan::ArgParseArgument::STRING, "METHOD");
	setDefaultValue(trimOpt, "WIN");
	setValidValues(trimOpt, "WIN BWA TAIL");
	addOption(parser, trimOpt);

	// ADAPTER TRIMMING
	addSection(parser, "Adapter removal options");
	seqan::ArgParseOption adapterFileOpt = seqan::ArgParseOption(
				"a", "adapters", "FastA file containing the two adapter sequences.",
				seqan::ArgParseArgument::STRING, "ADAPTER_FILE");
	addOption(parser, adapterFileOpt);

	seqan::ArgParseOption pairedModeOpt = seqan::ArgParseOption(
				"p", "paired", "Trim two input files with paired-end algorithm.",
				seqan::ArgParseOption::INTEGER, "BOOL");
	setDefaultValue(pairedModeOpt, 1);
	addOption(parser, pairedModeOpt);

	seqan::ArgParseOption rateOpt = seqan::ArgParseOption(
			"e", "errors", "Allowed errors in adapter detection. Double in [0,1) "
					       "is error rate, larger integers are absolute errors.",
			seqan::ArgParseOption::DOUBLE, "VALUE");
	addOption(parser, rateOpt);

	seqan::ArgParseOption overlapOpt = seqan::ArgParseOption(
				"o", "overlap", "Minimum length of overlap for a significant adapter match.",
				seqan::ArgParseOption::INTEGER, "VALUE");
		addOption(parser, overlapOpt);

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

	// Parse quality trimming parameters.
	seqan::CharString method;
	TrimmingAlgorithm trim;
	int cutoff, min_length;

	getOptionValue(method, parser, "m");
	if (method == "WIN")
		trim = Mean(5);
	else if (method == "BWA")
		trim = BWA();
	else
		trim = Tail();

	getOptionValue(cutoff, parser, "q");
	getOptionValue(min_length, parser, "l");

	// Parse adapter trimming parameters.
	bool paired = fileCount == 2;
	getOptionValue(paired, parser, "p");

	seqan::CharString adapterFile;
	seqan::Dna5String adapter1, adapter2, id;

	// If adapter sequences are given, we read them in any case.
	if (isSet(parser, "a"))
		getOptionValue(adapterFile, parser, "a");

		seqan::SequenceStream adapterStream(toCString(adapterFile));
		if (!seqan::isGood(adapterStream)){
			// TODO: Does this really mean "not found"?
			std::cerr << "Adapter file '" << adapterFile << "' not found. \n";
			return 1;
		}

		int err = seqan::readRecord(id, adapter1, adapterStream);
			err += seqan::readRecord(id, adapter2, adapterStream);
		if (err != 0){
			std::cerr << "Error while reading adapters from '" << adapterFile << "'.\n";
			return 1;
		}
	// If they are not given, but we would need them (single-end trimming), output error.
	else if (!paired)
	{
		std::cerr << "Unpaired adapter removal requires adapter sequences.\n";
		return 1;
	}

	Mode m;
	// User must fully specify mode, if he wants to. (But not specifying both is ok too.)
	if (seqan::isSet(parser, "e") != seqan::isSet(parser, "o")){
		std::cerr << "User must define both error rate and minimum overlap.\n";
		return 1;
	}

	// Read both options (if one is given, the above check guarantees that the other is too.)
	if (seqan::isSet(parser, "e")){
		int o;
		double e;
		getOptionValue(o, parser, "o");
		getOptionValue(e, parser, "e");
		m = User(o, round(e)); // TODO: Make behave as promised by option description or change description.
	}
	// Otherwise use the automatic configuration.
	else
		m = Auto();

	// TODO: We now have the quality and adapter trimming options. They are mostly mutually valid
	//       for both single and paired-end trimming. Now pass them into specific methods to deal
	//       with actual single- and paired-end processing.

	/*
	 seqan::CharString fileName1;
	 getArgumentValue(fileName1, parser, 0, 0);

	 if (fileCount == 1)
	   blabla for one file.
	 else
	   blabla for two files.

	// Main program loop for processing batches.
	 */


    return 0;
}

