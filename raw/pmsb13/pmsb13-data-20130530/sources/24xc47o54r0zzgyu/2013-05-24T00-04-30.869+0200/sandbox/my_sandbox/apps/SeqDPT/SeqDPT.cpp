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
#define SEQAN_PROFILE
#ifdef SEQAN_ENABLE_DEBUG
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

// Headers for creating subdirectories.
#include <errno.h>
// For setting the number of threads used by OpenMP.
#include <omp.h>

#include "readTrimming.h"
#include "adapterTrimming.h"
#include "demultiplex.h"

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef seqan::String<seqan::Dna5Q> Dna5QString;

/*!
 * \brief Enum for storing the trimming mode selected by the user.
 */
enum TrimmingMode
{
	E_WINDOW,
	E_BWA,
	E_TAIL
};

/*!
 * \brief Enum for storing the match options used for adapter trimming.
 */
enum MatchMode
{
	E_AUTO,
	E_USER
};

/*!
\brief Struct holding all demultiplexing parameters.
*/
struct DemultiplexingParams
{
	seqan::String<char> barcodeFile;						/*!< Holds the path to the barcode-file*/
	seqan::StringSet<seqan::String<seqan::Dna> > barcodes;	/*!< Holds the StringSet of barcodes*/
	seqan::StringSet<seqan::String<char> > barcodeIds;		/*!< Holds the StringSet of barcode-IDs*/
	seqan::String<char> multiplexFile;						/*!< Holds the path to the multiplex-file*/
	seqan::StringSet<seqan::String<seqan::Dna5Q> > multiplex; /*!< Holds the StringSet of multiplex barcodes*/
	bool approximate;										/*!< TRUE if approximate search shall be used*/
	bool hardClip;											/*!< TRUE if hardClip shall be used*/
	bool run;												/*!< TRUE if demultiplexing shall run*/
	bool runx;												/*!< TRUE if multiplex demultiplexing shall run*/

	DemultiplexStats stats;									/*!< Holds interesting numbers about demultiplexing.*/
};
/*!
\brief Struct holding all adapter trimming parameters.
*/
/*!
 * \brief Struct holding the parameters needed for adapter trimming.
 */
struct AdapterTrimmingParams
{
	bool paired;		/*!< TRUE if paired-end data is used.*/
	bool noAdapter;		/*!< TRUE if no adapter file is provided.*/
	bool run;			/*!< TRUE if the adapter trimming shall be run.*/
	Dna5QString adapter1;	/*!< Holds the first adapter.*/
	Dna5QString adapter2;	/*!< Holds the second adapter*/
	Mode mode;				/*!< Hold the ::Mode-object which determines which trimming mode shall be used.*/
	MatchMode mmode;        /*!< Enum, indicating which trimming mode is used. Necessary for casting mode. */

	AdapterTrimmingStats stats; /*!< Holds interesting numbers about trimmed adapters.*/

	AdapterTrimmingParams() : paired(false), noAdapter(false), run(false), mmode(E_AUTO) {};
};
/*!
\brief Struct holding all quality trimming  parameters.
*/
/*!
 * \brief Struct holding the parameters needed for quality trimming.
 */
struct QualityTrimmingParams
{
	TrimmingMode trim_mode;			/*!< Holds the ::TrimmingAlgorithm-object which determines which algrotihm shall be used.*/
	int cutoff;						/*!< Holds the cutoff score.*/
	int min_length;					/*!< Holds the minam length of a sequnce after trimming.*/
	bool run;						/*!< TRUE if the quality trimming shall run.*/

	QualityTrimmingStats stats; 	/*!< Holds interesting numbers about quality trimming.*/

	QualityTrimmingParams() : trim_mode(E_WINDOW), cutoff(-1), min_length(1), run(false) {};
};

/*!
 * \brief Struct that hold program parameters that might need to be
 * available at multiple places in the program.
 */
struct ProgramParams
{
	int fileCount;
	int readCount;
	double processTime, ioTime;
	seqan::SequenceStream fileStream1, fileStream2;

	ProgramParams() : fileCount(0), readCount(0), processTime(0), ioTime(0) {};
};

/*!
 * \brief Class that dynamically manages output streams that write out
 * sets of sequences.
 */
class OutputStreams
{
	typedef seqan::SequenceStream * PSeqStream;
	typedef seqan::Pair<PSeqStream, PSeqStream> TStreamPair;
	std::map<int, TStreamPair> pairedFileStreams;
	std::map<int, PSeqStream> fileStreams;

	const seqan::CharString basePath;
	seqan::CharString extension;

public:
	/*!
	 * \brief Constructor for the OutputStreams object. Prepares the file extension which will
	 * be used for all streams created by this object and saves a base directory path.
	 */
	OutputStreams(seqan::CharString base, seqan::SeqIOFileFormat_::Type format, bool compress) : basePath(base)
	{
		seqan::CharString fileExt("");
		switch (format)
		{
			case seqan::SeqIOFileFormat_::FILE_FORMAT_FASTQ:
				seqan::append(fileExt, seqan::CharString(".fastq"));
				break;
			case seqan::SeqIOFileFormat_::FILE_FORMAT_FASTA:
				seqan::append(fileExt, seqan::CharString(".fasta"));
				break;
			default:
			{
				std::cerr << "File format error. << std::endl";
				seqan::append(fileExt, seqan::CharString(".fasta"));
			}
		}

		if (compress)
			seqan::append(fileExt, seqan::CharString(".gz"));

		extension = fileExt;
	}

	/*!
	 * \brief Checks whether a key exists in a std::map.
	 * \return True if the key exists, false otherwise.
	 */
	template <typename TKey, typename TMap>
	bool exists(TKey& key, TMap& map)
	{
		return map.find(key) == map.end();
	}

	/*!
	 * \brief Add a new output streams to the collection of streams.
	 * \param fileName The name of the first file that will be created.
	 * \param id The associated id that will be used to identify the stream.
	 */
	void addStream(seqan::CharString fileName, int id)
	{
		seqan::CharString path = basePath;
		seqan::append(path, fileName);
		seqan::append(path, this->extension);
		char* file = seqan::toCString(path);

		PSeqStream stream = new SequenceStream(file, seqan::SequenceStream::WRITE);
		fileStreams[id] = stream;
	}

	/*!
	 * \brief Add a new output streams to the collection of streams.
	 * \param fileName1 The name of the first file that will be created.
	 * \param fileName2 The name of the second file that will be created.
	 * \param id The associated id that will be used to identify the pair of streams.
	 */
	void addStreams(seqan::CharString fileName1, seqan::CharString fileName2, int id)
	{
		// Prepend basePath and append file extension to the filename.
		seqan::CharString path1 = basePath, path2 = basePath;
		seqan::append(path1, fileName1);
		seqan::append(path1, this->extension);
		seqan::append(path2, fileName2);
		seqan::append(path2, this->extension);

		char* file1 = seqan::toCString(path1);
		char* file2 = seqan::toCString(path2);

		PSeqStream stream1 = new SequenceStream(file1, seqan::SequenceStream::WRITE);
		PSeqStream stream2 = new SequenceStream(file2, seqan::SequenceStream::WRITE);

		pairedFileStreams[id] = TStreamPair(stream1, stream2);
	}

	/*!
	 * \brief This method takes a vector of numbers and checks if these numbers are
	 * already associated with a stream. If not, a new stream is added and the opened
	 * file is named according to the list of names. One or two files are created.
	 * \param map The list of IDs for which the existence of a file should be checked.
	 * \param names The list of names to be used when creating new files.
	 * \param pair Indicates whether one or two (a pair of files) should be created.
	 */
	template <typename TMap, typename TNames>
	void updateStreams(TMap& map, TNames& names, bool pair)
	{
		for (unsigned i=0; i < length(map); ++i)
		{
			unsigned streamIndex = map[i];
			bool missing = pair ? exists(streamIndex, pairedFileStreams) : exists(streamIndex, fileStreams);
			// If no stream for this id exists, create one.
			if (missing)
			{
				// If the index is 0 (unidentified) create special stream.
				// Otherwise use index to get appropriate name for output file.
				seqan::CharString file;
				if (streamIndex > 0)
					file = names[streamIndex-1];
				else
					file = seqan::CharString("unidentified");

				// Add file extension to stream and create it.
				if (pair)
				{

// We use a linux syscall here. TODO: Try on Windows with _mkdir (should be compatible).
#ifdef __linux__
					// Create a new subfolder at basePath/[barcodeID, unidentified].
					seqan::CharString folderPath(basePath);
					seqan::append(folderPath, file);
					if (mkdir(seqan::toCString(folderPath), 0777) == -1)
					{
						if (errno == EEXIST)
							std::cerr << "Warning: Directory " << folderPath << " already exists.\n";
					}

					// Turn file target from [barcodeID,unidentified]
					// to subfolder [barcodeID, unidentified]/[barcodeID, unidentified]
					seqan::CharString filePath(file);
					seqan::append(file, "/");
					seqan::append(file, filePath);
#endif

					// To differentiate between the paired reads, add index to the file name.
					seqan::CharString file2 = file;
					seqan::append(file, seqan::CharString("_1"));
					seqan::append(file2, seqan::CharString("_2"));
					this->addStreams(file, file2, streamIndex);
				}
				else
					this->addStream(file, streamIndex);
			}
		}
	}

	/*!
	 * \brief Writes the sets of ids and sequences to their corresponding
	 * files.
	 * \param ids A list of sets of sequence IDs. (As returned by readRecord etc.)
	 * \param seqs A list of sets of sequences.
	 * \param map A mapping of the sets of sequences to their corresponding output streams.
	 * \param names Names to be used when creating new streams.
	 */
	template <typename TIds, typename TSeqs, typename TMap, typename TNames>
	void writeSeqs(TIds& ids, TSeqs& seqs, TMap& map, TNames& names)
	{
		updateStreams(map, names, false);
		for (unsigned i=0; i < length(seqs); ++i)
		{
			unsigned streamIndex = map[i];
			seqan::writeAll(*fileStreams[streamIndex], ids[i], seqs[i]);
		}
	}

	/*!
	 * \brief Writes the sets of ids and sequences to their corresponding
	 * files. Overload for writing paired-end sequence sets.
	 * \param ids1 A list of sets of forward sequence IDs. (As returned by readRecord etc.)
	 * \param seqs1 A list of sets of forward sequences.
	 * \param ids2 A list of sets of backward sequence IDs. (As returned by readRecord etc.)
	 * \param seqs2 A list of sets of backward sequences.
	 * \param map A mapping of the sets of sequences to their corresponding output streams.
	 * \param names Names to be used when creating new streams.
	 */
	template <typename TIds, typename TSeqs, typename TMap, typename TNames>
	void writeSeqs(TIds& ids1, TSeqs& seqs1, TIds& ids2, TSeqs& seqs2, TMap& map, TNames& names)
	{
		updateStreams(map, names, true);
		for (unsigned i=0; i < length(seqs1); ++i)
		{
			unsigned streamIndex = map[i];
			TStreamPair tmp = pairedFileStreams[streamIndex];
			seqan::writeAll(*tmp.i1, ids1[i], seqs1[i]);
			seqan::writeAll(*tmp.i2, ids2[i], seqs2[i]);
		}
	}

	 /*!
	  * \brief Destructor of the object holding the output streams. Needed
	  * to destroy the streams properly after they have been created with new.
	  */
	~OutputStreams()
	{
		// Delete all created file streams.
		for (unsigned i=0; i < length(fileStreams); ++i)
			delete fileStreams[i];

		for (unsigned i=0; i < length(pairedFileStreams); ++i)
		{
			TStreamPair tmp = pairedFileStreams[i];
			delete tmp.i1;
			delete tmp.i2;
		}
	}
};

// ============================================================================
// Functions
// ============================================================================

/*!
\brief Function for loading the sequence-files.
\param seqStream SequenceStream-object of the file.
\param ids StringSet of String of chars the IDs shall be stored in.
\param seqs StringSet of Dna5Q-Strings the sequences shall be stored in.
\param records Unsigned int determining the number of records to be read in one go.
\param params The ::DemultiplexingParams-object.
\return An integer: \b 1 on errors, \b 0 otherwise.
*/
int loadSeqs(seqan::SequenceStream& seqStream, seqan::StringSet<seqan::String<char> >& ids, seqan::StringSet<seqan::String<seqan::Dna5Q> >& seqs, unsigned records)
{
	resize(ids, records, Exact());
	resize(seqs, records, Exact());
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
/*!
\brief Function for loading all barcodes.
\param path Pointer to the barcode file.
\param params The ::DemultiplexingParams-object the barcodes shall be written to.
\return An interger: \b 1 on errors, \b otherwise.
*/
int loadBarcodes(char const * path, DemultiplexingParams& params)
{
	seqan::SequenceStream bcStream(path, seqan::SequenceStream::READ);

	if(!seqan::isGood (bcStream))
	{
		std::cerr << "Error while opening file'" <<  params.barcodeFile << "'.\n";
		return 1;
	}

	if(seqan::readAll(params.barcodeIds, params.barcodes, bcStream) != 0)
	{
		std::cerr << "Error while reading barcodes from '" << params.barcodeFile << "'.\n";
		return 1;
	}
	resize(params.stats.groups, length(params.barcodes));	//Sets the right size for the stats vector and fills it with 0.		
	for(unsigned i = 0; i < length(params.stats.groups); ++i)
		params.stats.groups[i] = 0;
	return 0;
}

// Kann glaube ich gel?scht werden...
int loadMultiplex (seqan::SequenceStream& multiplexStream, DemultiplexingParams& params, unsigned records)
{
	resize(params.multiplex, records, Exact());
	seqan::StringSet<seqan::String<char> > ids;
	resize(ids, records);
	if (!seqan::isGood(multiplexStream))
	{
		std::cerr << "Error while opening file '" <<  params.multiplexFile << "'.\n";
		return 1;
	}
	if(seqan::readBatch(ids, params.multiplex, multiplexStream, records) != 0)
	{
		std::cerr << "Error while reading barcodes from '" << params.multiplexFile << "'.\n";
		return 1;
	}
	return 0;
}

int openStream(seqan::CharString const & file, seqan::SequenceStream & stream)
{
	open(stream, seqan::toCString(file));
	if (!isGood(stream))
	{
		std::cerr << "Error while opening input file '" << file << "'.\n";
		return 1;
	}

	return 0;
}

/*!
\brief Function for loading all demultiplexing parameters.
\param parser The SeqAn parser-object.
\param params The ::DemultiplexingParams-object the paramters shall be written to.
\return An interger: \b 1 on errors, \b 0 otherwise.
*/
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

	// APPROXIMATE/EXACT MATCHING---------------------
	params.approximate = seqan::isSet(parser, "app");
	// HARD CLIP MODE --------------------------------
	params.hardClip = seqan::isSet(parser, "hc");
	// RUN
	params.run = isSet(parser, "b");
	if (isSet(parser, "x"))
		params.runx = true;
	else
		params.runx = false;
	return 0;
}
/*!
\brief Function for loading all adapter trimming parameters.
\param parser The SeqAn parser-object.
\param params The ::AdapterTrimmingParams-object parameters shall be written to.
\return An interger: \b 1 on errors, \b 0 otherwise.
*/
int loadAdapterTrimmingParams(seqan::ArgumentParser const& parser, AdapterTrimmingParams & params)
{
	// PAIRED-END ------------------------------
	int fileCount =  getArgumentValueCount(parser, 0);
	// Only consider paired-end mode if two files are given and user wants paired mode.
	params.paired = fileCount == 2 && !isSet(parser, "np");
	params.noAdapter = isSet(parser, "na");

	// Set run flag, depending on essential parameters.
	params.run = isSet(parser,"a") || (params.noAdapter && fileCount == 2);

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
	else if (isSet(parser, "np") || (params.noAdapter && fileCount == 1))
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
		if (!params.run)
		{
			std::cerr << "Adapter removal alignment parameters require adapters or --no-adapter flag.\n";
			return 1;
		}

		int o;
		int e;
		getOptionValue(o, parser, "o");
		getOptionValue(e, parser, "e");
		params.mode = User(o, e);
		params.mmode = E_USER;
	}
	// Otherwise use the automatic configuration.
	else
	{
		params.mode = Auto();
		params.mmode = E_AUTO;
	}

	return 0;
}
/*!
\brief Function for loading all quality trimming parameters.
\param parser The SeqAn parser-object.
\param params The ::QualityTrimmingParams-object the parameters shall be written to.
\return An interger: \b 1 on errors, \b 0 otherwise.
*/
int loadQualityTrimmingParams(seqan::ArgumentParser const & parser, QualityTrimmingParams & params)
{
	// TRIMMING METHOD ----------------------------
	std::string method;
	getOptionValue(method, parser, "m");
	if (method == "WIN")
		params.trim_mode = E_WINDOW;
	else if (method == "BWA")
		params.trim_mode = E_BWA;
	else
		params.trim_mode = E_TAIL;
	// QUALITY CUTOFF ----------------------------
	if (isSet(parser, "q"))
		getOptionValue(params.cutoff, parser, "q");
	// MINIMUM SEQUENCE LENGTH -------------------
	getOptionValue(params.min_length, parser, "l");

	// Set run flag, depending on essential parameters. (Which are in a valid state at this point.)
	params.run = isSet(parser, "q");
	return 0;
}

/*!
 * \brief Loads file names and
 */
int loadProgramParams(seqan::ArgumentParser const & parser, ProgramParams & params)
{
	params.fileCount = getArgumentValueCount(parser, 0);

	// Load files.
	seqan::CharString fileName1, fileName2;

	getArgumentValue(fileName1, parser, 0, 0);
	if (openStream(fileName1, params.fileStream1) != 0)
		return 1;

	if (params.fileCount == 2)
	{
		getArgumentValue(fileName2, parser, 0, 1);
		if (openStream(fileName2, params.fileStream2) != 0)
			return 1;

		if (params.fileStream1._fileFormat != params.fileStream2._fileFormat)
		{
			std::cerr << "Input files should have the same file format.\n";
			return 1;
		}
	}

	return 0;
}

int checkParams(ProgramParams const & programParams, DemultiplexingParams const & demultiplexingParams,
		AdapterTrimmingParams const & adapterTrimmingParams, QualityTrimmingParams & qualityTrimmingParams)
{
	// Were there options that activated at least one processing stage?
	if (!(adapterTrimmingParams.run || qualityTrimmingParams.run || demultiplexingParams.run))
	{
		std::cerr << "No processing stage was specified.\n";
		return 1;
	}

	// If quality trimming was selected, check if file format includes qualities.
	if (qualityTrimmingParams.run)
	{
		if (programParams.fileStream1._fileFormat != seqan::SeqIOFileFormat_::FILE_FORMAT_FASTQ
		  || ((programParams.fileCount == 2) && programParams.fileStream2._fileFormat != seqan::SeqIOFileFormat_::FILE_FORMAT_FASTQ))
		{
			std::cerr << "Quality trimming requires quality information, please specify fastq files." << std::endl;
			return 1;
		}
	}

	// If we don't demultiplex (and therefore take file names from the barcode IDs), set file names.
	if (!demultiplexingParams.run)
	{

		//std::cout << basename(toCString(fileName1)) << "\n\n";
	}

	return 0;
}

// PROGRAM STAGES ---------------------
/*!
\brief Function for checking the parameters and calling the desired demultiplexing functions
\param params The ::DemultiplexingParams object.
\param seqs StringSet of sequences the operations shall be performed on.
\param ids StringSet of IDs of the sequences.
\param esaFinder the Finder-object holding information on the barcod-index.
\param map Map the information on the barcode groups shall be stored in.
\return An integer: \b 0 on errors, \b 1 otherwise.
*/
// DEMULTIPLEXING

template <typename TSeqsVec, typename TIdsVec, typename TFinder,typename TMulti, typename TMap>
int demultiplexingStage(DemultiplexingParams& params, TSeqsVec& seqs, TIdsVec& ids, TFinder& esaFinder, TMulti& multiplexStream, TMap& map, unsigned records)
{
	if (!params.run)
		return 0;

	std::vector<std::vector<int> > groups;
	if (params.runx && !params.approximate)
	{
		//DEBUG_MSG(std::cout << "Demultiplexing exact multiplex single-end reads.\n");
		if(loadMultiplex(multiplexStream, params, records)!=0)
			return 1;
		groups = DoAll(params.multiplex, params.barcodes, esaFinder, params.stats);
	}
	else if(!params.approximate)
	{
		if(!check(seqs[0], ids[0], params.barcodes))		// On Errors with barcodes return 1;
			return 1;					
		//DEBUG_MSG("Demultiplexing exact inline single-end reads.\n");
		groups = DoAll(seqs[0], params.barcodes, esaFinder, params.hardClip, params.stats);
	}
	else if(params.runx && params.approximate)
	{
		//DEBUG_MSG("Demultiplexing approximate multiplex single-end reads.\n");
		if(loadMultiplex(multiplexStream, params, records)!=0)
			return 1;
		groups = DoAll(params.multiplex, params.barcodes, params.stats);
	}
	else
	{
		if(!check(seqs[0], ids[0], params.barcodes))			// On Errors with barcodes return 1;
			return 1;
		//DEBUG_MSG("Demultiplexing approximate inline single-end reads.\n");
		groups = DoAll(seqs[0], params.barcodes, params.hardClip, params.stats);
	}

	// Saves information on how groups correspond to barcodes.
	map.clear();
	for (unsigned i = 0; i < length(groups); ++i)
		if (length(groups[i]) != 0) map.push_back(i);

	// Sorting the results into the sequence- and ID vectors
	std::vector<seqan::StringSet<seqan::String<seqan::Dna5Q> > > sortedSeqs; 
	std::vector<seqan::StringSet<seqan::String<char> > > sortedIds;
	buildSets(seqs[0], ids[0], groups, sortedSeqs, sortedIds);

	resize(seqs, length(sortedSeqs));
	resize(ids, length(sortedIds));
	seqs = sortedSeqs;
	ids = sortedIds;
	return 0;
}
/*!
\brief Overload of ::demultiplexingStage(const DemultiplexingParams& params, TSeqsVec& seqs, TIdsVec& ids, TFinder& esaFinder, TMap& map)
\param params The ::DemultiplexingParams object.
\param seqs StringSet of forward-reads the operations shall be performed on.
\param seqsRev StringSet of backward-reads the operations shall be performes on.
\param ids StringSet of IDs of the sequences.
\param esaFinder the Finder-object holding information on the barcod-index.
\param map Map the information on the barcode groups shall be stored in.
\return An integer: \b 0 on errors, \b 1 otherwise.
*/

template <typename TSeqsVec, typename TIdsVec, typename TFinder, typename TMulti, typename TMap>
int demultiplexingStage(DemultiplexingParams& params, TSeqsVec& seqs, TSeqsVec& seqsRev, TIdsVec& ids, TIdsVec& idsRev, TFinder& esaFinder, TMulti& multiplexStream, TMap& map, unsigned records)
{
	if (!params.run)
		return 0;

	std::vector<std::vector<int> > groups;
	if (params.runx && !params.approximate)
	{
		//DEBUG_MSG("Demultiplexing exact multiplex paired-end reads.\n");
		if(loadMultiplex(multiplexStream, params, records)!=0)
			return 1;
		groups = DoAll(params.multiplex, params.barcodes, esaFinder, params.stats);
	}
	else if(!params.approximate)
	{
		if(!check(seqs[0], seqsRev[0],  ids[0], params.barcodes))		// On Erros with barcodes return 1;
			return 1;					
		//DEBUG_MSG("Demultiplexing exact inline paired-end reads.\n");
		groups = DoAll(seqs[0], params.barcodes, esaFinder, params.hardClip, params.stats);
	}
	else if(params.runx && params.approximate)
	{
		//DEBUG_MSG("Demultiplexing approximate multiplex paired-end reads.\n");
		if(loadMultiplex(multiplexStream, params, records)!=0)
			return 1;
		groups = DoAll(params.multiplex, params.barcodes, params.stats);
	}
	else
	{
		if(!check(seqs[0], seqsRev[0], ids[0], params.barcodes))			// On Erros with barcodes return 1;
			return 1;
		//DEBUG_MSG("Demultiplexing approximate inline paired-end reads.\n");
		groups = DoAll(seqs[0], params.barcodes, params.hardClip, params.stats);
	}

	// Saves information on how groups correspond to barcodes.
	map.clear();
	for (unsigned i = 0; i < length(groups); ++i)
		if (length(groups[i]) != 0) map.push_back(i);

	// Sorting the results into the sequence- and ID vectors
	std::vector<seqan::StringSet<seqan::String<seqan::Dna5Q> > > sortedSeqs; 
	std::vector<seqan::StringSet<seqan::String<seqan::Dna5Q> > > sortedSeqsRev; 
	std::vector<seqan::StringSet<seqan::String<char> > > sortedIds;
	std::vector<seqan::StringSet<seqan::String<char> > > sortedIdsRev;

	buildSets(seqs[0], seqsRev[0] ,ids[0], idsRev[0], groups, sortedSeqs, sortedSeqsRev, sortedIds, sortedIdsRev);

	resize(seqs, length(sortedSeqs));
	resize(seqsRev, length(sortedSeqs));
	resize(ids, length(sortedIds));
	resize(idsRev, length(sortedIds));

	seqs = sortedSeqs;
	seqsRev = sortedSeqsRev;
	ids = sortedIds;
	idsRev = sortedIdsRev;
	return 0;
}

// ADAPTER TRIMMING

/*!
 * \brief Executes adapter trimming functions for single-end file input.
 * \param params Adapter trimming parameters to be used for trimming.
 * \param seqSet The set of sequences to be trimmed.
 */
template <typename TSeqs>
void adapterTrimmingStage(AdapterTrimmingParams& params, TSeqs& seqSet)
{
	if (!params.run)
		return;

	//DEBUG_MSG("Trimming single-end adapters.\n");

	typename seqan::Iterator<TSeqs, Rooted>::Type it;

	// Templates don't support runtime polymorphism, so code paths for all possibilities.
	switch(params.mmode)
	{
	case E_USER:
		for (it = seqan::begin(seqSet); it != seqan::end(seqSet); it++)
				stripAdapterBatch(value(it), params.adapter2, (User&) params.mode, params.stats);
		break;
	case E_AUTO:
		for (it = seqan::begin(seqSet); it != seqan::end(seqSet); it++)
				stripAdapterBatch(value(it), params.adapter2, (Auto&) params.mode, params.stats);
		break;
	}
}

/*!
 * \brief Executes adapter trimming functions for paired-end file input.
 * \param params Adapter trimming parameters to be used for trimming.
 * \param seqSet1 The set of forward sequences to be trimmed.
 * \param seqSet2 The set of backward sequences to be trimmed.
 */
template <typename TSeqs>
void adapterTrimmingStage(AdapterTrimmingParams& params, TSeqs& seqSet1, TSeqs& seqSet2)
{
	if (!params.run)
		return;

	typename seqan::Iterator<TSeqs, Rooted>::Type it1, it2;
	for (it1 = seqan::begin(seqSet1), it2 = seqan::begin(seqSet2); it1 != seqan::end(seqSet1); it1++, it2++)
	{
		if (!params.paired)
		{
			//DEBUG_MSG("Trimming paired-end adapters in single-end mode.\n");
			// Templates don't support runtime polymorphism, so code paths for all possibilities.
			switch(params.mmode)
			{
			case E_USER:
			{
				stripAdapterBatch(value(it1), params.adapter2, (User&)params.mode, params.stats);
				stripReverseAdapterBatch(value(it2), params.adapter1, (User&)params.mode, params.stats);
				break;
			}
			case E_AUTO:
			{
				stripAdapterBatch(value(it1), params.adapter2, (Auto&)params.mode, params.stats);
				stripReverseAdapterBatch(value(it2), params.adapter1, (Auto&)params.mode, params.stats);
				break;
			}
			}
		}
		else
		{
			//DEBUG_MSG("Trimming paired-end adapters.\n");
			stripPairBatch(value(it1), value(it2), params.stats);
		}
	}
}

// QUALITY TRIMMING

/*!
 * \brief Executes quality trimming functions for single-end file input.
 * \param params Quality trimming parameters to be used for trimming.
 * \param idSet The ids of the set of sequences to be trimmed.
 * \param seqSet The set of sequences to be trimmed.
 * \remark Including the IDs even though only the sequences are trimmed is necessary
 * in case we remove some sequence. We also have to remove the corresponding ID in that case.
 */
template <typename TIds, typename TSeqs>
void qualityTrimmingStage(QualityTrimmingParams& params, TIds& idSet, TSeqs& seqSet)
{
	if (params.run)
	{
		//DEBUG_MSG("Trimming qualities.\n");
		// Templates don't support runtime polymorphism, so code paths for all possibilities.
		switch (params.trim_mode)
		{
		case E_WINDOW:
			for (unsigned i=0; i < length(idSet); ++i)
				trimBatch(seqSet[i], params.cutoff, Mean(5));
			break;
		case E_BWA:
			for (unsigned i=0; i < length(idSet); ++i)
				trimBatch(seqSet[i], params.cutoff, BWA());
			break;
		case E_TAIL:
			for (unsigned i=0; i < length(idSet); ++i)
				trimBatch(seqSet[i], params.cutoff, Tail());
			break;
		}
	}

	for (unsigned i=0; i < length(idSet); ++i)
		dropReads(idSet[i], seqSet[i], params.min_length, params.stats);
}

/*!
 * \brief Executes quality trimming functions for paired-end file input.
 * \param params Quality trimming parameters to be used for trimming.
 * \param idSet1 The ids of the set of forward sequences to be trimmed.
 * \param seqSet1 The set of forward sequences to be trimmed.
 * \param idSet2 The ids of the set of backward sequences to be trimmed.
 * \param seqSet2 The set of backward sequences to be trimmed.
 * \remark Including the IDs even though only the sequences are trimmed is necessary
 * in case we remove some sequence. We also have to remove the corresponding ID in that case.
 */
template <typename TIds, typename TSeqs>
void qualityTrimmingStage(QualityTrimmingParams& params, TIds& idSet1, TSeqs& seqSet1, TIds& idSet2, TSeqs& seqSet2)
{
	if (params.run)
	{
		//DEBUG_MSG("Trimming (pair) qualities.\n");
		// Templates don't support runtime polymorphism, so code paths for all possibilities.
		switch (params.trim_mode)
		{
		case E_WINDOW:
			for (unsigned i=0; i < length(idSet1); ++i)
				trimPairBatch(seqSet1[i], seqSet2[i], params.cutoff, Mean(5));
			break;
		case E_BWA:
			for (unsigned i=0; i < length(idSet1); ++i)
				trimPairBatch(seqSet1[i], seqSet2[i], params.cutoff, BWA());
			break;
		case E_TAIL:
			for (unsigned i=0; i < length(idSet1); ++i)
				trimPairBatch(seqSet1[i], seqSet2[i], params.cutoff, Tail());
			break;
		}
	}

	for (unsigned i=0; i < length(idSet1); ++i)
		dropReads(idSet1[i], seqSet1[i], idSet2[i], seqSet2[i], params.min_length, params.stats);

}

// END PROGRAM STAGES ---------------------

void printStatistics(ProgramParams& programParams, DemultiplexingParams& demultiplexParams, 
				AdapterTrimmingParams& adapterParams, QualityTrimmingParams& qualityParams)
{
	bool paired = programParams.fileCount == 2;
	bool adapter = adapterParams.run;
	//bool quality = qualityParams.run;

	int read_factor = (1+paired);

	std::cout << std::endl;
	std::cout << "Read statistics\n";
	std::cout << "===============\n";

	std::cout << "Processed reads: " << read_factor*programParams.readCount;
	if (paired) std::cout << " (2 * " << programParams.readCount << ")";
	std::cout << std::endl;
	std::cout << "  Dropped reads: " << qualityParams.stats.dropped_1 + qualityParams.stats.dropped_2 << "\n\n";

	//Statistics for Demultiplexing
	unsigned usedBarcodes = 0;
	if(demultiplexParams.run)	
	{
		std::cout << "Barcode Demultiplexing statistics\n";
		std::cout << "=================================\n";
		
		for (unsigned i = 0; i < length(demultiplexParams.stats.groups); ++i)
		{
			if(demultiplexParams.stats.groups[i] != 0)
				++usedBarcodes;
		}
		if(usedBarcodes > 0)	//correction for unidentified group
			--usedBarcodes;
	}

	std::cout << "Barcodes used: " << usedBarcodes << " of " << length(demultiplexParams.barcodes) << "\n";
	std::cout << "Reads per Barcode:\n";
	std::cout << "Unidentified: " << demultiplexParams.stats.groups[0] << "(" << (double)demultiplexParams.stats.groups[0]/(double)programParams.readCount*100 << "%)\n";
	for (unsigned i = 1; i < length(demultiplexParams.stats.groups); ++i)
	{
		std::cout << demultiplexParams.barcodeIds[i]<<": " << demultiplexParams.stats.groups[i] << "(" << (double)demultiplexParams.stats.groups[i]/(double)programParams.readCount*100 << "%)\n";
	}
	std::cout << std::endl;

	std::cout << "File statistics\n";
	std::cout << "===============\n";

	// How many reads are left.
	int survived1 = programParams.readCount - qualityParams.stats.dropped_1;
	int survived2 = programParams.readCount - qualityParams.stats.dropped_2;

	// In percentage points.
	double surv_proc_1 = (double)survived1/(double)programParams.readCount*100;
	double surv_proc_2 = (double)survived2/(double)programParams.readCount*100;

	std::cout << "File 1:\n";
	std::cout << "-------\n";
	std::cout << "  Surviving: " << survived1 << "/" << programParams.readCount
			  << " (" << std::setprecision(3) << surv_proc_1 << "%)\n";
	if (adapter) std::cout << "   Adapters: " << adapterParams.stats.a1count << "\n";
	std::cout << std::endl;

	if (paired)
	{
		std::cout << "File 2:\n";
		std::cout << "-------\n";
		std::cout << "  Surviving: " << survived2 << "/"
				  << programParams.readCount << " (" << surv_proc_2 << "%)\n";
		if (adapter) std::cout << "   Adapters: " << adapterParams.stats.a2count << "\n";
		std::cout << std::endl;
	}


	if (adapter)
	{
		int mean = adapterParams.stats.overlapSum/(read_factor*programParams.readCount);
		std::cout << "Adapter sizes:\n";
		std::cout << "Min: " << adapterParams.stats.minOverlap << ",  Mean: " << mean
				<< ", Max: " << adapterParams.stats.maxOverlap << "\n\n";
	}

	// Print processing and IO time. IO is (approx.) the whole loop without the processing part.
	std::cout << "Timming statistics:\n";
	std::cout << "==================\n";
	std::cout << "Processing time: " << std::setw(5) << programParams.processTime << " seconds.\n";
	std::cout << "       I/O time: " << std::setw(5) << programParams.ioTime << " seconds.\n";
	std::cout << std::endl;
}

// END FUNCTION DEFINITIONS ---------------------------------------------


/*!
 * \brief This method defines the argument parser for SeqDPT.
 */
seqan::ArgumentParser initParser(){
	// PARSER DEFINITON ----------------------------------------
	seqan::ArgumentParser parser("SeqDPT");
	setShortDescription(parser, "The Sequence Processing Toolkit");
	addUsageLine(parser, " READ_FILE1 [READ_FILE2] [OPTIONS]");
	addDescription(parser,
            "SeqDPT is a tool for processing of sequenced NGS reads. It "
			"is possible to demultiplex the reads and order them according to different kind of barcodes, to remove adapter "
			"contamination from reads and to trim low quality bases. The different tools are controlled through "
			"command line parameters and can operate on both single- and paired-end read data.");
	setDate(parser, __DATE__);
	setVersion(parser, "1.95.2.4.1 (SeqDPT 2013 - Academic license)");
	seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUTFILE, "FILES", true);
	setValidValues(fileArg, "fasta fasta.gz fastq fastq.gz");
	addArgument(parser, fileArg);

	// GENERAL OPTIONS -----------------------------------------
	addSection(parser, "General Options");
	
	seqan::ArgParseOption recordOpt = seqan::ArgParseOption(
		"r", "records", "Number of records to be read in one run.",
		seqan::ArgParseOption::INTEGER, "VALUE");
	setDefaultValue(recordOpt, 1000);
	addOption(parser, recordOpt);

	seqan::ArgParseOption compressOpt = seqan::ArgParseOption(
			"c", "compress", "Compress output files with gzip.");
		addOption(parser, compressOpt);

	seqan::ArgParseOption threadOpt = seqan::ArgParseOption(
				"tnum", "threads", "Number of threads used.",
				seqan::ArgParseOption::INTEGER, "THREADS");
	setDefaultValue(threadOpt, 1);
	setMinValue(threadOpt, "1");
	addOption(parser, threadOpt);

	seqan::ArgParseOption outputOpt = seqan::ArgParseOption(
		"out", "output", "Folder for output (must already exist).",
		seqan::ArgParseOption::STRING, "OUTPUT");
	setDefaultValue(outputOpt, "");
	addOption(parser, outputOpt);

	// Barcode Demultiplexing
	addSection(parser, "Demultiplexing Options");
	
	seqan::ArgParseOption barcodeFileOpt = seqan::ArgParseOption(
		"b", "barcodes", "FastA file containing the used barcodes and their IDs. Necessary for demutiplexing.",
		seqan::ArgParseArgument::STRING, "BARCODE_FILE");
	addOption(parser, barcodeFileOpt);

	seqan::ArgParseOption multiplexFileOpt = seqan::ArgParseOption(
		"x", "multiplex", "FastA/FastQ file containing the barcode for each read.",
		seqan::ArgParseArgument::STRING, "MULTIPLEX_FILE");
	addOption(parser, multiplexFileOpt);
	
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
				"na", "no-adapter", "Trim adapters from paired-end reads without using reference adapters.");
	addOption(parser, noAdapterOpt);

	seqan::ArgParseOption pairedModeOpt = seqan::ArgParseOption(
				"np", "no-paired", "Trim paired-end input with single-end trimming method.");
	addOption(parser, pairedModeOpt);

	seqan::ArgParseOption rateOpt = seqan::ArgParseOption(
			"e", "errors", "Allowed errors in adapter detection.",
			seqan::ArgParseOption::INTEGER, "VALUE");
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

	seqan::CharString output;
	getOptionValue(output, parser, "out");

	int threads;
	getOptionValue(threads, parser, "tnum");
	omp_set_num_threads(threads);

	//--------------------------------------------------
	// Parse demultiplexing parameters.
	//--------------------------------------------------

	DemultiplexingParams demultiplexingParams;

	getOptionValue(demultiplexingParams.multiplexFile, parser, "x");
	seqan::SequenceStream multiplexStream;	//Initialising the SequenceStream for the multiplex file
	if (isSet(parser, "x"))
		open(multiplexStream, seqan::toCString(demultiplexingParams.multiplexFile));
			
	if(loadDemultiplexingParams(parser, demultiplexingParams) != 0)
		return 1;

	//--------------------------------------------------
	// Process Barcodes
	//--------------------------------------------------

	seqan::Index<seqan::StringSet<seqan::String<seqan::Dna> >, seqan::IndexEsa<> > indexSet(demultiplexingParams.barcodes);
	seqan::Finder<seqan::Index<seqan::StringSet<seqan::String<seqan::Dna> >, seqan::IndexEsa<> > > esaFinder(indexSet);

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
	// General program parameters and additional checks.
	//--------------------------------------------------

	ProgramParams programParams;
	if ( loadProgramParams(parser, programParams) != 0)
		return 1;

	if (checkParams(programParams, demultiplexingParams, adapterTrimmingParams, qualityTrimmingParams) != 0)
		return 1;

	//--------------------------------------------------
	// Processing
	//--------------------------------------------------

	// Prepare output stream object and initial mapping from StringSets to files.
	OutputStreams outputStreams(output, programParams.fileStream1._fileFormat, isSet(parser, "c"));
	std::vector<unsigned> map(1,0);

	// Start processing. Different functions are needed for one or two input files.
	std::cout << "\nProcessing reads." << std::endl;

	SEQAN_PROTIMESTART(loopTime);
	if (fileCount == 1)
	{
		if (!demultiplexingParams.run)
			outputStreams.addStream(seqan::CharString("processedReads"), 0);

		std::vector<seqan::StringSet<seqan::CharString> > idSet;
		std::vector<seqan::StringSet<Dna5QString> > seqSet;

		while (!atEnd(programParams.fileStream1))
		{
			idSet.clear();
			seqSet.clear();

			idSet.push_back(seqan::StringSet<seqan::CharString>());
			seqSet.push_back(seqan::StringSet<Dna5QString>());

			if (loadSeqs(programParams.fileStream1, idSet[0], seqSet[0], records) == 0)
			{
				programParams.readCount += length(idSet[0]);
				SEQAN_PROTIMESTART(processTime);

				// Demultiplexing
				if(demultiplexingStage(demultiplexingParams, seqSet, idSet,
									esaFinder, multiplexStream, map, records) != 0)
					return 1;

				// Adapter trimming.
				adapterTrimmingStage(adapterTrimmingParams, seqSet);

				// Quality trimming.
				qualityTrimmingStage(qualityTrimmingParams, idSet, seqSet);
				programParams.processTime += SEQAN_PROTIMEDIFF(processTime);

				// Append to output file.
				outputStreams.writeSeqs(idSet, seqSet, map, demultiplexingParams.barcodeIds);

			}
			else
				return 1;
		}
	 }
	 else
	 {
		if (!demultiplexingParams.run)
			outputStreams.addStreams(seqan::CharString("processedReads_1"), seqan::CharString("processedReads_2"), 0);

		std::vector<seqan::StringSet<seqan::CharString> > idSet1, idSet2;
		std::vector<seqan::StringSet<Dna5QString> > seqSet1, seqSet2;

		while (!(atEnd(programParams.fileStream1) || atEnd(programParams.fileStream2)))
		{
			idSet1.clear(); idSet2.clear();
			seqSet1.clear(); seqSet2.clear();

			idSet1.push_back(seqan::StringSet<seqan::CharString>());
			seqSet1.push_back(seqan::StringSet<Dna5QString>());
			idSet2.push_back(seqan::StringSet<seqan::CharString>());
			seqSet2.push_back(seqan::StringSet<Dna5QString>());

			if (loadSeqs(programParams.fileStream1, idSet1[0], seqSet1[0], records) == 0 &&
				loadSeqs(programParams.fileStream2, idSet2[0], seqSet2[0], records) == 0)
			{
				programParams.readCount += length(idSet1[0]);
				SEQAN_PROTIMESTART(processTime); // Measure processing time.

				// Demultiplexing.
				if(demultiplexingStage(demultiplexingParams, seqSet1, seqSet2, idSet1, idSet2,
									esaFinder, multiplexStream, map, records) != 0)
					return 1;

				// Adapter trimming.
				adapterTrimmingStage(adapterTrimmingParams, seqSet1, seqSet2);

				// Quality trimming.
				qualityTrimmingStage(qualityTrimmingParams, idSet1, seqSet1, idSet2, seqSet2);
				programParams.processTime += SEQAN_PROTIMEDIFF(processTime); // End of processing time.

				// Append to output file.
				outputStreams.writeSeqs(idSet1, seqSet1, idSet2, seqSet2,
										map, demultiplexingParams.barcodeIds);
			}
			else return 1;
		}
	}

	double loop = SEQAN_PROTIMEDIFF(loopTime);
	programParams.ioTime = loop - programParams.processTime;

	printStatistics(programParams, demultiplexingParams, adapterTrimmingParams, qualityTrimmingParams);

	return 0;
}
