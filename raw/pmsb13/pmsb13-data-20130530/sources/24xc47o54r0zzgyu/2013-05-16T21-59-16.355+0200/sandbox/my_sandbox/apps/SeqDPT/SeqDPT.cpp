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

typedef seqan::String<seqan::Dna5Q> Dna5QString;

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

	seqan::ArgParseOption recordOpt = seqan::ArgParseOption(
		"r", "records", "Number of records to be read in one run.",
		seqan::ArgParseOption::INTEGER, "VALUE");
	setDefaultValue(recordOpt, 1000);
	addOption(parser,recordOpt);

	seqan::ArgParseOption outputOpt = seqan::ArgParseOption(
		"out", "output", "Folder for output (must allready exist)",
		seqan::ArgParseOption::STRING, "OUTPUT");
	addOption(parser,outputOpt);

	// Barcode Demultiplexing
	addSection(parser, "Demultiplexing Options");
	
	seqan::ArgParseOption barcodeFileOpt = seqan::ArgParseOption(
		"b", "barcodes", "FastA file containing the used barcodes and their IDs. Necessary for demutiplexing.",
		seqan::ArgParseArgument::STRING, "BARCODE_FILE");
	addOption(parser,barcodeFileOpt);

	seqan::ArgParseOption multiplexFileOpt = seqan::ArgParseOption(
		"x", "multiplex", "FastA/FastQ file containing the barcode for each read.",
		seqan::ArgParseArgument::STRING, "MULTIPLEX_FILE");
	addOption(parser,multiplexFileOpt);
	
	addOption(parser, seqan::ArgParseOption(
        "app", "approximate", "Select approximate barcode demultiplexing, allowing one mismatch."));

	addOption(parser, seqan::ArgParseOption(
        "hc", "hardClip", "Select hardClip option for barcode clipping, clipping the first length(barcode) bases even if no matching barcode has been found."));

	// READ TRIMMING
	addSection(parser, "Quality trimming options");
	seqan::ArgParseOption qualOpt = seqan::ArgParseOption(
			"q", "quality", "Quality threshold for read trimming.",
			seqan::ArgParseArgument::INTEGER, "PHRED");
	//setRequired(qualOpt, true);
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

	// Parse general paramters.
	unsigned records = 1000;
	getOptionValue(records, parser, "r");
	seqan::String<char> output;
	getOptionValue(output, parser, "out");

	// Parse demultiplexing parameters.
	bool approximate = isSet(parser, "approximate");
	bool hardClip = isSet(parser, "hardClip");
	seqan::CharString barcodeFile;
	seqan::StringSet<seqan::String<char> >  barcodeIds;
	seqan::StringSet<seqan::String<seqan::Dna5Q> > barcodes;
	seqan::CharString multiplexFile;
	seqan::StringSet<seqan::String<seqan::Dna5Q> > multiplex;
	if (isSet(parser, "b"))
	{
		getOptionValue(barcodeFile, parser, "b");
		seqan::SequenceStream barcodeStream(toCString(barcodeFile), SequenceStream::READ);

		if (!seqan::isGood(barcodeStream)){
			std::cerr << "Error while opening the barcode-file.\n" ;
			return 1;
		}

		int err = seqan::readAll(barcodeIds, barcodes, barcodeStream);
		if (err != 0){
			std::cerr << "Error while reading the barcodes.\n";
			return 1;
		}
		if (isSet(parser, "x"))
		{
			getOptionValue(multiplexFile, parser, "x");
			seqan::SequenceStream multiplexStream(toCString(multiplexFile), SequenceStream::READ);
													
			StringSet<String<char> > muids;
			StringSet<String<char> > muquals;
			resize(multiplex, records);
			resize(muids, records);
			resize(muquals, records);

			if (!isGood(multiplexStream))
			{
				std::cerr << "Error while opening the multiplex barcode-file.\n";
				return 1;
			}
			if(readBatch(muids, multiplex, muquals, multiplexStream, records) != 0)		
			{
				std::cerr << "Error while reading the multiplex barcodes.\n";
				return 1;
			}
			return 0;
		}
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
	// Only consider paired-end mode if two files are given.
	if (paired)
		getOptionValue(paired, parser, "p");

	seqan::CharString adapterFile, id;
	Dna5QString adapter1, adapter2;

	// If adapter sequences are given, we read them in any case.
	if (isSet(parser, "a"))
	{
		getOptionValue(adapterFile, parser, "a");

		seqan::SequenceStream adapterStream(toCString(adapterFile));
		if (!seqan::isGood(adapterStream)){
			std::cerr << "Error while opening file'" << adapterFile << "'.\n";
			return 1;
		}

		int err = seqan::readRecord(id, adapter1, adapterStream);
			err += seqan::readRecord(id, adapter2, adapterStream);
		if (err != 0){
			std::cerr << "Error while reading adapters from '" << adapterFile << "'.\n";
			return 1;
		}
	}
	// If they are not given, but we would need them (single-end trimming), output error.
	else if (!paired)
	{
//		std::cerr << "Unpaired adapter removal requires adapter sequences.\n";
//		return 1;
	}

	Mode mode;
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
		mode = User(o, round(e)); // TODO: Make behave as promised by option description or change description.
	}
	// Otherwise use the automatic configuration.
	else
		mode = Auto();

	// TODO: We now have the quality and adapter trimming options. They are mostly mutually valid
	//       for both single and paired-end trimming. Now pass them into specific methods to deal
	//       with actual single- and paired-end processing.


	seqan::CharString fileName1;
	getArgumentValue(fileName1, parser, 0, 0);
	seqan::SequenceStream fileStream1(toCString(fileName1));
	if (!isGood(fileStream1))
	{
		std::cerr << "Error while opening input file " << fileName1 << ".\n";
		return 1;
	}

	seqan::SequenceStream outStream1("test_output1.fq", seqan::SequenceStream::WRITE);

	if (fileCount == 1)
	{
		seqan::StringSet<seqan::CharString> ids;
		seqan::StringSet<Dna5QString> seqs;
		seqan::StringSet<seqan::String<char> > quals;

		while (!atEnd(fileStream1))
		{
			if (seqan::readBatch(ids, seqs, fileStream1, records) == 0)
			{
				// Demultiplexing
				if(isSet(parser, "b"))
				{
					if(isSet(parser, "x"))
					{
						if(isSet(parser, "q"))			//TODO Das "a" muss evtl noch durch einen anderen paramter ersetzt werden.
						{
							if(isSet(parser, "a"))		//Multiplex Demultiplexing with Adapter Removal and Read Trimming
							{
								//....
							}
							else						//Multiplex Demultiplexing with Read Trimming
							{
								//...					
							}
						}
						else							
						{
							if(isSet(parser, "a"))		//Multiplex Demultiplexing with Adapter Removal
							{
								//....
							}
							else						//Multiplex Demultiplexing
							{
								std::cout << "Test.\n";
								std::vector<std::vector<int> > groups = DoAll(seqs, ids, multiplex, barcodes, approximate);
								std::vector<StringSet<String<Dna5Q> > > gSeqs;
								std::vector<StringSet<String<char> > > gIds;
								buildSets(seqs, ids, groups, gSeqs, gIds);
								writeGroups(gSeqs, gIds, barcodeIds, groups, output);
							}
						}
					}
					else
					{
						if(isSet(parser, "q"))
						{
							if(isSet(parser, "a"))		//Inline Demultiplexing with Adapter Removal and Read Trimming
							{
								//...
							}
							else						//Inline Demultiplexing with Read Trimming
							{
								//...
							}
						}
						else
						{
							if(isSet(parser, "a"))		//Inline Demultiplexing with Adapter Removal
							{
								//...
							}
							else						//Inline Demultiplexing
							{
								std::vector<std::vector<int> > groups = DoAll(seqs, ids, barcodeFile, approximate, hardClip);
								std::vector<StringSet<String<Dna5Q> > > gSeqs;
								std::vector<StringSet<String<char> > > gIds;
								buildSets(seqs, ids, groups, gSeqs, gIds);
								writeGroups(gSeqs, gIds, barcodeIds, groups, output);
							}
						}
					}
				}

				// Adapter trimming.

		//		stripAdapterBatch(seqs, adapter2, Auto());

				// Quality trimming.
				//Mean m = Mean(5);
		//		trimBatch(ids, seqs, cutoff, min_length, Mean(5));
				//std::cout << trimmed << "\n";
				// Append to output file.
		//		seqan::writeAll(outStream1, ids, seqs);
			}
			else
			{
				std::cerr << "Error while reading input file " << fileName1 << ".\n";
				return 1;
			}
		}
	 }
	 else
	 {
		seqan::CharString fileName2;
		getArgumentValue(fileName2, parser, 0, 1);

		seqan::SequenceStream fileStream2(toCString(fileName2));
		if (!isGood(fileStream2))
		{
			std::cerr << "Error while opening input file " << fileName2 << ".\n";
			return 1;
		}

		seqan::StringSet<seqan::CharString> ids1, ids2;
		seqan::StringSet<Dna5QString> seqs1, seqs2;

		while (!(atEnd(fileStream1) || atEnd(fileStream2)))
		{
			if (seqan::readBatch(ids1, seqs1, fileStream1, records) == 0 &&
				seqan::readBatch(ids2, seqs2, fileStream1, records) == 0)
			{
				// Demultiplexing.

				// Adapter trimming.

				// Quality trimming.

				// Append to output file.
			}
			else
			{
				std::cerr << "Error while reading input files\n";
				return 1;
			}
		}
	}
	return 0;
}

