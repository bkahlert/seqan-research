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


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef seqan::String<seqan::Dna5Q> Dna5QString;

/*!
 Struct holding all demultiplexing paramters.
*/
struct DemultiplexingParams
{
	seqan::String<char> barcodeFile;						/*! Holds the path to the barcode-file*/
	seqan::StringSet<seqan::String<seqan::Dna> > barcodes;	/*! Holds the StringSet of barcodes*/
	seqan::StringSet<seqan::String<seqan::Dna> > barcodeIds;/*! Holds the StringSet of barcode-IDs*/
	seqan::String<char> multiplexFile;						/*! Holds the path to the multiplex-file*/
	seqan::StringSet<seqan::String<seqan::Dna> > multiplex; /*! Holds the StringSet of multiplex barcodes*/
	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet;	/*! Holds the IndexSet for exact search*/
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder; /*! Holds the Finder for above indexSet*/
	bool approximate;										/*! TRUE if approximate search shall be used*/
	bool hardClip;											/*! TRUE if hardClip shall be used*/
	bool run;												/*! TRUE if demultiplexing shall run*/
	bool runx;												/*! TRUE if multiplex demultiplexing shall run*/
};

struct AdapterTrimmingParams
{
	bool paired, noAdapter, run;
	Dna5QString adapter1, adapter2;
	Mode mode;

	AdapterTrimmingParams() : paired(false), noAdapter(false), run(false) {};
};

struct QualityTrimmingParams
{
	TrimmingAlgorithm trim_method;
	int cutoff;
	int min_length;
	bool run;

	QualityTrimmingParams() : cutoff(-1), min_length(1), run(false) {};
};

struct OutputStreams
{
	typedef seqan::SequenceStream * PSeqStream;
	std::vector<seqan::Pair<PSeqStream, PSeqStream> > fileStreams;
};

// FUNCTION DEFINITIONS -------------------------------------------------

template<typename TParams>
int loadDemultiplexingParams(seqan::ArgumentParser const& parser, TParams& params)
{
	// BARCODES -----------------------------------
	if (isSet(parser, "b"))
	{
		getOptionValue(params.barcodeFile, parser, "b");
		if(loadBarcodes(toCString(params.barcodeFile), params) != 0) 
			return 1;
	}
	// APPROXIMATE MATCHING---------------------------
	params.approximate = seqan::isSet(parser, "app");
	if(!params.approximate)
	{	//Creating Index and Finder for barcodes.
		params.indexSet = Index<StringSet<String<Dna> >, IndexEsa<> > params.indexSet(params.barcodes);
		params.esaFinder = Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > params.esaFinder(params.indexSet);
	}
	// HARD CLIP MODE --------------------------------
	params.hardClip = seqan::isSet(parser, "hc");
	// RUN
	params.run = true;
	return 0;
}
//Overload for multiplex barcoding
template<typename TParams, typename TMultiplexStream>
int loadDemultiplexingParams(seqan::ArgumentParser const& parser, TParams& params, TMultiplexStream& multiplexStream, unsigned records)
{
	// BARCODES, APPROXIMATE, HARDCLIP--------------
	if (loadDemultiplexingParams(parser, params) != 0)
		return 1;

	// MULTIPLEX FILE -----------------------------
	getOptionValue(params.multiplexFile, parser, "x");
	if(loadMultiplex(multiplexStream, params, records) != 0) 
		return 1;
	params.runx = true;
	return 0;
}

int loadAdapterTrimmingParams(seqan::ArgumentParser const& parser, AdapterTrimmingParams & params)
{
	// PAIRED-END ------------------------------
	int fileCount =  getArgumentValueCount(parser, 0);
	// Only consider paired-end mode if two files are given and user wants paired mode.
	params.paired = fileCount == 2 && !isSet(parser, "np");
	params.noAdapter = isSet(parser, "na");

	// ADAPTER SEQUENCES ----------------------------
	seqan::CharString adapterFile, id;

	// If adapter sequences are given, we read them in any case.
	if (isSet(parser, "a"))
	{
		getOptionValue(adapterFile, parser, "a");

		seqan::SequenceStream adapterStream(toCString(adapterFile));
		if (!seqan::isGood(adapterStream)){
			std::cerr << "Error while opening file'" << adapterFile << "'.\n";
			return 1;
		}

		int err = seqan::readRecord(id, params.adapter1, adapterStream);
			err += seqan::readRecord(id, params.adapter2, adapterStream);
		if (err != 0){
			std::cerr << "Error while reading adapters from '" << adapterFile << "'.\n";
			return 1;
		}
	}
	// If they are not given, but we would need them (single-end trimming), output error.
	else if (!params.paired)
	{
		std::cerr << "Unpaired adapter removal requires adapter sequences.\n";
		return 1;
	}

	// TRIMMING MODE ----------------------------
	// User must fully specify mode, if he wants to. (But not specifying both is ok too.)
	if (seqan::isSet(parser, "e") != seqan::isSet(parser, "o")){
		std::cerr << "User must define both error rate (-e) and minimum overlap (-o).\n";
		return 1;
	}

	// Read both options (if one is given, the above check guarantees that the other is too.)
	if (seqan::isSet(parser, "e")){
		// If user tried to specify alignment parameters, but didn't activate adapter trimming, warn and exit.
		if (!params.noAdapter && !isSet(parser, "a"))
		{
			std::cerr << "Adapter removal alignment parameters require adapters or --no-adapter flag.\n";
			return 1;
		}

		int o;
		double e;
		getOptionValue(o, parser, "o");
		getOptionValue(e, parser, "e");
		params.mode = User(o, round(e)); // TODO: Make behave as promised by option description or change description.
	}
	// Otherwise use the automatic configuration.
	else
		params.mode = Auto();

	// Set run flag, depending on essential parameters. (Which are in a valid state at this point.)
	params.run = isSet(parser,"a") || isSet(parser, "na");
	return 0;
}

int loadQualityTrimmingParams(seqan::ArgumentParser const & parser, QualityTrimmingParams & params)
{
	// TRIMMING METHOD ----------------------------
	std::string method;
	getOptionValue(method, parser, "m");
	if (method == "WIN")
		params.trim_method = Mean(5);
	else if (method == "BWA")
		params.trim_method = BWA();
	else
		params.trim_method = Tail();
	// QUALITY CUTOFF ----------------------------
	if (isSet(parser, "q"))
		getOptionValue(params.cutoff, parser, "q");
	// MINIMUM SEQUENCE LENGTH -------------------
	getOptionValue(params.min_length, parser, "l");

	// Set run flag, depending on essential parameters. (Which are in a valid state at this point.)
	params.run = isSet(parser, "q");
	return 0;
}

int loadSeqs(seqan::SequenceStream seqStream, seqan::StringSet<seqan::String<char> >& ids, seqan::StringSet<seqan::String<seqan::Dna5Q> >& seqs, unsigned records)
{
	resize(ids, records);
	resize(seqs, records);
	if (!seqan::isGood(seqStream))
	{
		std::cerr << "Error while opening the sequence-file.\n";
		return 1;
	}
	
	if(seqan::readBatch(ids, seqs, seqStream, records) != 0)		
	{
		std::cerr << "Error while reading the sequences.\n";
		return 1;
	}
	return 0;
}

int loadBarcodes(char const * path, DemultiplexingParams& demultiplexingParams)
{
	seqan::SequenceStream bcStream(path, seqan::SequenceStream::READ);
	
	if(!seqan::isGood (bcStream))
	{
		std::cerr << "Error while opening file'" <<  demultiplexingParams.barcodeFile << "'.\n";
		return 1;
	}

	if(seqan::readAll(demultiplexingParams.barcodeIds, demultiplexingParams.barcodes, bcStream) != 0)
	{
		std::cerr << "Error while reading barcodes from '" << demultiplexingParams.barcodeFile << "'.\n";
		return 1;
	}
	return 0;
}

int loadMultiplex (seqan::SequenceStream multiplexStream, DemultiplexingParams& demultiplexingParams, unsigned records)
{
	resize(demultiplexingParams.multiplex, records);
	seqan::StringSet<seqan::String<char> > ids;
	resize(ids, records);
	if (!seqan::isGood(multiplexStream))
	{
		std::cerr << "Error while opening file'" <<  demultiplexingParams.multiplexFile << "'.\n";
		return 1;
	}
	if(seqan::readBatch(ids, demultiplexingParams.multiplex, multiplexStream, records) != 0)		
	{
		std::cerr << "Error while reading barcodes from '" << demultiplexingParams.multiplexFile << "'.\n";
		return 1;
	}
	return 0;
}

// PROGRAM STAGES

template <typename TSeqsVec, typename TIdsVec>
int DemultiplexStage(const DemultiplexingParams& params, TSeqsVec& seqs, TIdsVec& Ids)
{
	if (!params.run)
		return 0;

	if (params.runx && !params.approximate)
	{
		std::cout << "Demultiplexing exact multiplex single-end reads.\n";
		std::vector<std::vector<int> > groups = DoAll(params.multiplex, param.barcodes, params.esaFinder);
	}
	else if(!params.approximate)
	{
		if(!check(seqs[0], ids[0], params.barcodes))		// On Erros with barcodes return 1;
			return 1;					
		std::cout << "Demultiplexing exact inline single-end reads.\n";
		std::vector<std::vector<int> > groups = DoAll(seqs[0], ids[0], params.barcodes, params.esaFinder, params.hardClip);
	}
	else if(params.runx && params.approximate)
	{
		std::cout << "Demultiplexing approximate multiplex single-end reads.\n";
		std::vector<std::vector<int> > groups = DoAll(params.multiplex, params.barcodes);
	}
	else
	{
		if(!check(seqs[0], ids[0], params.barcodes))			// On Erros with barcodes return 1;
			return 1;
		std::cout << "Demultiplexing approximate inline single-end reads.\n";
		std::vector<std::vector<int> > groups = DoAll(seqs[0], ids[0], params.barcodes, params.hardClip);
	}
	// Sorting the results into the sequence vector
	std::vector<TSeqs>& sortedSeqs; 
	std::vector<TIds>& sortedIds;
	buildSets(seqs[0], ids[0], groups, sortedSeqs, sortedIds);
	resize(seqs, length(sortedSeqs));
	resize(ids, length(sortedIds));
	for(unsigned i = 0; i < length(sortedSeqs); ++i)
	{
		for(unsigned j = 0; j < length(sortedSeqs[i]); ++j)
		{
			appendValue(seqs[i], seqs[sortedSeqs[i][j]]);
			appendValue(ids[i],ids[sortedIds[i][j]]);
		}
	}
	return 0;
}

template <typename TSeqs>
void adapterTrimmingStage(AdapterTrimmingParams& params, TSeqs& seqSet)
{
	stripAdapterBatch(seqSet, params.adapter2, Auto());
}

template <typename TSeqs>
void adapterTrimmingStage(AdapterTrimmingParams& params, TSeqs& seqSet1, TSeqs& seqSet2)
{
	if (!params.run)
		return;

	if (!params.paired)
	{
		std::cout << "Trimming paired-end adapters in single-end mode.\n";
		stripAdapterBatch(seqSet1, params.adapter2, Auto());
		stripReverseAdapterBatch(seqSet2, params.adapter1, Auto());
	}
	else
	{
		std::cout << "Trimming paired-end adapters.\n";
		stripPairBatch(seqSet1, seqSet2);
	}
}

template <typename TIds, typename TSeqs>
void qualityTrimmingStage(QualityTrimmingParams& params, TIds& idSet, TSeqs& seqSet)
{
	if (!params.run)
		return;

	std::cout << "Trimming qualities.\n";
	trimBatch(idSet, seqSet, params.cutoff, params.min_length, Mean(5));
}

template <typename TIds, typename TSeqs>
void qualityTrimmingStage(QualityTrimmingParams& params, TIds& idSet1, TSeqs& seqSet1, TIds& idSet2, TSeqs& seqSet2)
{
	if (!params.run)
		return;

	std::cout << "Trimming qualities.\n";
	trimPairBatch(idSet1, seqSet1, idSet2, seqSet2, params.cutoff, params.min_length, Mean(5));
}

// END FUNCTION DEFINITIONS ---------------------------------------------




seqan::ArgumentParser initParser(){
	// PARSER DEFINITON ----------------------------------------
	seqan::ArgumentParser parser("SeqDPT");
	addDescription(parser,
                   "Program for inline/multiplex barcode demultipexing, adapter removal and low-quality read trimming. For single-end and paired-end reads.");
	setDate(parser, __DATE__);
	addArgument(parser, seqan::ArgParseArgument(
			seqan::ArgParseArgument::INPUTFILE, "FILES", true));

	// GENERAL OPTIONS -----------------------------------------
	addSection(parser, "General Options");
	
	seqan::ArgParseOption recordOpt = seqan::ArgParseOption(
		"r", "records", "Number of records to be read in one run.",
		seqan::ArgParseOption::INTEGER, "VALUE");
	setDefaultValue(recordOpt, 1000);
	addOption(parser,recordOpt);

	seqan::ArgParseOption outputOpt = seqan::ArgParseOption(
		"out", "output", "Folder for output (must allready exist).",
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

	seqan::ArgParseOption noAdapterOpt = seqan::ArgParseOption(
				"na", "no-adapter", "Trim adapters from paired-end reads without using adapter sequences.");
	addOption(parser, noAdapterOpt);

	seqan::ArgParseOption pairedModeOpt = seqan::ArgParseOption(
				"np", "no-paired", "Trim paired-end input with single-end trimming method.");
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

	return parser;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------
// Program entry point.

int main(int argc, char const ** argv)
{
	seqan::ArgumentParser parser = initParser();

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

	//--------------------------------------------------
	// Parse general parameters.
	//--------------------------------------------------

	unsigned records;
	getOptionValue(records, parser, "r");
	seqan::String<char> output;
	getOptionValue(output, parser, "out");

	//--------------------------------------------------
	// Parse demultiplexing parameters.
	//--------------------------------------------------

	DemultiplexingParams demultiplexingParams;
	if (isSet(parser, "x"))
	{
		seqan::SequenceStream multiplexStream(toCString(demultiplexingParams.multiplexFile), seqan::SequenceStream::READ);
		if(loadDemultiplexingParams(parser, demultiplexingParams, multiplexStream, records) != 0)
			return 1;
	}
	else if(loadDemultiplexingParams(parser, demultiplexingParams) != 0)
		return 1;

	//--------------------------------------------------
	// Parse quality trimming parameters.
	//--------------------------------------------------

	QualityTrimmingParams qualityTrimmingParams;
	if ( loadQualityTrimmingParams(parser, qualityTrimmingParams) != 0 )
		return 1;

	//--------------------------------------------------
	// Parse adapter trimming parameters.
	//--------------------------------------------------

	AdapterTrimmingParams adapterTrimmingParams;
	if ( loadAdapterTrimmingParams(parser, adapterTrimmingParams) != 0 )
		return 1;

	//--------------------------------------------------
	// Processing
	//--------------------------------------------------

	seqan::CharString fileName1;
	getArgumentValue(fileName1, parser, 0, 0);
	seqan::SequenceStream fileStream1(toCString(fileName1));
	if (!isGood(fileStream1))
	{
		std::cerr << "Error while opening input file " << fileName1 << ".\n";
		return 1;
	}

	seqan::SequenceStream outStream1("test_output1.fq", seqan::SequenceStream::WRITE);

	if (!(adapterTrimmingParams.run || qualityTrimmingParams.run))
	{
		std::cerr << "No processing stage was specified.\n";
		return 1;
	}

	if (fileCount == 1)
	{
		seqan::StringSet<seqan::CharString> ids;
		seqan::StringSet<Dna5QString> seqs;

		while (!atEnd(fileStream1))
		{
			if (seqan::readBatch(ids, seqs, fileStream1, records) == 0)
			{
				

				// Adapter trimming.
				adapterTrimmingStage(adapterTrimmingParams, seqs);

				// Quality trimming.
				qualityTrimmingStage(qualityTrimmingParams, ids, seqs);

				// Append to output file.
				seqan::writeAll(outStream1, ids, seqs);
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

		seqan::SequenceStream outStream2("test_output2.fq", seqan::SequenceStream::WRITE);

		seqan::StringSet<seqan::CharString> ids1, ids2;
		seqan::StringSet<Dna5QString> seqs1, seqs2;

		while (!(atEnd(fileStream1) || atEnd(fileStream2)))
		{
			if (seqan::readBatch(ids1, seqs1, fileStream1, records) == 0 &&
				seqan::readBatch(ids2, seqs2, fileStream2, records) == 0)
			{
				// Demultiplexing.

				// Adapter trimming.
				adapterTrimmingStage(adapterTrimmingParams, seqs1, seqs2);
				// Quality trimming.
				qualityTrimmingStage(qualityTrimmingParams, ids1, seqs1, ids2, seqs2);

				// Append to output file.
				seqan::writeAll(outStream1, ids1, seqs1);
				seqan::writeAll(outStream2, ids2, seqs2);
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

