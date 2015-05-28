/*
 *   Options.h
 *
 *   Created on: Jul 31, 2012
 *       Author: jtr
 */

#ifndef FLEXBAR_OPTIONS_H_
#define FLEXBAR_OPTIONS_H_

#include <string>
#include <iostream>

#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>
#include <seqan/arg_parse.h>

#include "Enums.h"
#include "AdapterLoader.h"


const std::string getFlexbarBanner(const std::string version){
	
	std::string banner = "";
	
	banner += "               ________          __              \n";
	banner += "              / ____/ /__  _  __/ /_  ____ ______\n";
	banner += "             / /_  / / _ \\| |/ / __ \\/ __ `/ ___/\n";
	banner += "            / __/ / /  __/>  </ /_/ / /_/ / /    \n";
	banner += "           /_/   /_/\\___/_/|_/_.___/\\__,_/_/     \n\n";
	
	banner += "Flexbar - flexible barcode and adapter removal, version ";
	banner += version;
	banner += "\nBioinformatics in Quantitative Biology at BIMSB, GPLv3\n";
	
	return banner;
}


const std::string getFlexbarCitation(){
	return "\nMatthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich: Flexbar — flexible barcode and adapter processing for next-generation sequencing platforms. MDPI Biology 2012, 1(3):895-905.\n";
}


void defineOptionsAndHelp(seqan::ArgumentParser &parser, const std::string version, const std::string date){
	
	using namespace seqan;
	
	setVersion(parser, version);
	setDate(parser, date);
	
	setShortDescription(parser, "flexible barcode and adapter removal");
	
	addUsageLine(parser, "\\fB-t\\fP target \\fB-f\\fP format \\fB-r\\fP reads [\\fB-b\\fP barcodes] [\\fB-a\\fP adapters] [options]");
	
	// addOption(parser, ArgParseOption("v", "version", "Displays program version."));
	addOption(parser, ArgParseOption("c", "cite", "Shows citation information."));
	
	addSection(parser, "Basic options");
	addOption(parser, ArgParseOption("n", "threads", "Number of threads.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("t", "target", "Prefix for output file names.", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("r", "reads", "Input file with reads, that may contain barcodes.", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("p", "reads2", "Second input file for paired read scenario.", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("f", "format", "Input format of reads: csfasta, csfastq, fasta, fastq, fastq-sanger, fastq-solexa, fastq-i1.3, fastq-i1.5, fastq-i1.8 (illumina 1.8+).", ArgParseArgument::STRING));
	
	addSection(parser, "Barcode detection");
	addOption(parser, ArgParseOption("b",  "barcodes", "Fasta file with barcodes, specify (br) to use separate barcode reads.", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("br", "barcode-reads", "Fasta or fastq file with separate barcode reads, if not within reads.", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("be", "barcode-trim-end", "Type of barcoding within reads, see section trim-end modes.", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("bn", "barcode-tail-length", "Number of bases for tail trim-end modes. Default: barcode length.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("bo", "barcode-min-overlap", "Minimum overlap of barcode and read. Default: barcode length.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("bt", "barcode-threshold", "Allowed mismatches and indels per 10 bases for barcode.", ArgParseArgument::DOUBLE));
	addOption(parser, ArgParseOption("bk", "barcode-keep", "Keep barcodes within reads instead of removal."));
	addOption(parser, ArgParseOption("bm", "barcode-match", "Alignment match score.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("bi", "barcode-mismatch", "Alignment mismatch score.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("bg", "barcode-gap", "Alignment gap score.", ArgParseArgument::INTEGER));
	
	addSection(parser, "Adapter removal");
	addOption(parser, ArgParseOption("a",  "adapters", "Fasta file with adapters, or barcodes to remove within reads.", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("as", "adapter-seq", "Single adapter sequence as alternative to adapters option.", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("ae", "adapter-trim-end", "Type of adapter removal, see section trim-end modes.", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("an", "adapter-tail-length", "Number of bases for tail trim-end modes. Default: adapter length.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("ao", "adapter-min-overlap", "Minimum overlap of adapter and read.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("at", "adapter-threshold", "Allowed mismatches and indels per 10 bases for adapter.", ArgParseArgument::DOUBLE));
	addOption(parser, ArgParseOption("am", "adapter-match", "Alignment match score.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("ai", "adapter-mismatch", "Alignment mismatch score.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("ag", "adapter-gap", "Alignment gap score.", ArgParseArgument::INTEGER));
	
	addSection(parser, "Filtering and trimming");
	addOption(parser, ArgParseOption("u", "max-uncalled", "Allowed uncalled bases (N or .) in reads.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("x", "pre-trim-left", "Trim specified number of bases on 5' end of reads before detection.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("y", "pre-trim-right", "Trim specified number of bases on 3' end of reads before detection.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("q", "pre-trim-phred", "Trim 3' end until specified or higher quality reached.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("z", "post-trim-length", "Trim to specified read length from 3' end after removal.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("m", "min-readlength", "Minimum read length to remain after removal.", ArgParseArgument::INTEGER));
	
	addSection(parser, "Logging and tagging");
	addOption(parser, ArgParseOption("l", "log-level", "Print valid sequence alignments of reads.", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("s", "single-reads", "Output single reads for partially too short paired reads."));
	addOption(parser, ArgParseOption("o", "fasta-output", "Prefer non-quality formats fasta and csfasta for output."));
	addOption(parser, ArgParseOption("i", "length-dist", "Write length distribution for read output files."));
	addOption(parser, ArgParseOption("e", "number-tags", "Replace read tags by ascending number to save space."));
	addOption(parser, ArgParseOption("g", "removal-tags", "Tag reads that are subject to adapter or barcode removal."));
	addOption(parser, ArgParseOption("d", "random-tags", "Random read tags at barcode or adapter positions with N."));
	
	addSection(parser, "Trim-end modes");
	addText(parser._toolDoc, "ANY: longer part of read remains", false);
	addText(parser._toolDoc, "LEFT: align <= read end, right part remains", false);
	addText(parser._toolDoc, "RIGHT: align >= read start, left part remains", false);
	addText(parser._toolDoc, "LEFT_TAIL: consider first n bases of reads in alignment", false);
	addText(parser._toolDoc, "RIGHT_TAIL: use only last n bases, see tail-length options", false);
	
	
	// hideOption(parser, "removal-tags");
	
	setDefaultValue(parser, "threads",        "1");
	setDefaultValue(parser, "max-uncalled",   "0");
	setDefaultValue(parser, "min-readlength", "18");
	
	setDefaultValue(parser, "barcode-trim-end",  "ANY");
	setDefaultValue(parser, "barcode-threshold", "2.0");
	setDefaultValue(parser, "barcode-match",     "1");
	setDefaultValue(parser, "barcode-mismatch", "-1");
	setDefaultValue(parser, "barcode-gap",      "-7");
	
	setDefaultValue(parser, "adapter-trim-end",    "RIGHT");
	setDefaultValue(parser, "adapter-min-overlap", "1");
	setDefaultValue(parser, "adapter-threshold",   "3.0");
	setDefaultValue(parser, "adapter-match",       "1");
	setDefaultValue(parser, "adapter-mismatch",   "-1");
	setDefaultValue(parser, "adapter-gap",        "-7");
	
	setValidValues(parser, "log-level", "ALL MOD TAB");
	
	addTextSection(parser, "EXAMPLES");
	addText(parser._toolDoc, "\\fBflexbar\\fP \\fB-t\\fP target \\fB-f\\fP fastq-i1.3 \\fB-r\\fP reads.fastq   \\fB-b\\fP bar.fasta \\fB-a\\fP adap.fasta", false);
	addText(parser._toolDoc, "\\fBflexbar\\fP \\fB-t\\fP target \\fB-f\\fP csfastq    \\fB-r\\fP reads.csfastq \\fB-a\\fP adapters.fasta \\fB-ae\\fP LEFT");
}


struct Options{
	
	std::string readsFile, readsFile2, targetName, format;
	std::string adapterFile, barcodeFile, barReadsFile, adapterSeq, log_level;
	std::string a_trim_end, b_trim_end;
	
	bool isColorSpace, useAdapterFile, useNumberTag, useRemovalTag;
	bool switch2Fasta, writeSingleReads, writeLengthDist, randTag;
	
	int cutLen_begin, cutLen_end, phred_preQual, cutLen_read, a_tail_len, b_tail_len;
	int maxUncalled, min_readLen, a_min_overlap, b_min_overlap, nThreads;
	int match, mismatch, gapCost, b_match, b_mismatch, b_gapCost;
	
	float a_threshold, b_threshold;
	
	flexbar::TrimEnd         end, b_end;
	flexbar::FileFormat      fformat;
	flexbar::QualityType     qual;
	flexbar::RunType         runType;
	flexbar::LogLevel        logLevel;
	flexbar::BarcodeDetect   barDetect;
	flexbar::AdapterRemoval  adapRm;
	
	tbb::concurrent_vector<TAdapter> adapters, barcodes;
	
	Options(){
		
		isColorSpace     = false;
		useAdapterFile   = false;
		useNumberTag     = false;
		useRemovalTag    = false;
		writeSingleReads = false;
		writeLengthDist  = false;
		switch2Fasta     = false;
		randTag          = false;
		
		cutLen_begin  = 0;
		cutLen_end    = 0;
		cutLen_read   = 0;
		phred_preQual = 0;
		a_tail_len    = 0;
		b_tail_len    = 0;
		b_min_overlap = 0;
		
		logLevel  = flexbar::NONE;
		barDetect = flexbar::BOFF;
		adapRm    = flexbar::AOFF;
    }
};


void printLocalTime(){
	time_t t_current;
	time(&t_current);
	printf("Local time:            %s\n", asctime(localtime(&t_current)));
}


void parseCommandLine(seqan::ArgumentParser &parser, std::string version, int argc, char const ** argv){
	
	using namespace std;
	
	using seqan::ArgumentParser;
	
	cout << endl;
	
	ArgumentParser::ParseResult res = parse(parser, argc, argv);
	
	if(res != ArgumentParser::PARSE_OK){
		
		if(isSet(parser, "help")){
			cout << "\nFurther documentation on: flexbar.sourceforge.net" << endl;
		}
		cout << endl;
		
		exit(res == ArgumentParser::PARSE_ERROR);
	}
	
	if(isSet(parser, "cite")){
		cout << getFlexbarBanner(version) << getFlexbarCitation() << endl;
		exit(0);
	}
	
	if(! isSet(parser, "target") || ! isSet(parser, "format") || ! isSet(parser, "reads")){
		printShortHelp(parser);
		cout << "\nPlease set required parameters.\n" << endl;
		exit(1);
	}
	
	cout << getFlexbarBanner(version) << "\n" << endl;
	printLocalTime();
}


void loadProgramOptions(Options &o, seqan::ArgumentParser &parser){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::ArgumentParser;
	
	
	// basic options
	
	getOptionValue(o.targetName, parser, "target");
	cout << "Target name:           " << o.targetName << endl;
	
	getOptionValue(o.format, parser, "format");
	cout << "File format:           " << o.format;
	
	if(o.format == "fasta"){
		o.fformat = FASTA;
		o.qual    = SANGER;
	}
	else if(o.format == "fastq"){
		o.fformat = FASTQ;
		o.qual    = SANGER;
		
		if(isSet(parser, "pre-trim-phred")){
			cerr << "\nSpecify precise fastq format for quality based trimming.\n" << endl;
			exit(1);
		}
	}
	else if(o.format == "fastq-sanger"){
		o.fformat = FASTQ;
		o.qual    = SANGER;
	}
	else if(o.format == "fastq-solexa"){
		o.fformat = FASTQ;
		o.qual    = SOLEXA;
	}
	else if(o.format == "fastq-i1.3" || o.format == "fastq-illumina1.3"){
		o.fformat = FASTQ;
		o.qual    = ILLUMINA13;
	}
	else if(o.format == "fastq-i1.5" || o.format == "fastq-illumina1.5"){
		o.fformat = FASTQ;
		o.qual    = ILLUMINA13;
	}
	else if(o.format == "fastq-i1.8" || o.format == "fastq-illumina1.8"){
		o.fformat = FASTQ;
		o.qual    = SANGER;
	}
	else if(o.format == "csfasta"){
		o.fformat = CSFASTA;
		o.qual    = SANGER;
		o.isColorSpace = true;
	}
	else if(o.format == "csfastq"){
		o.fformat = CSFASTQ;
		o.qual    = SANGER;
		cout << " (sanger quality scaling)";
		o.isColorSpace = true;
	}
	else{
		cerr << "\n\n" << "Specified input file format is unknown!" << "\n" << endl;
		exit(1);
	}
	cout << endl;
	
	
	if(isSet(parser, "fasta-output") && (o.fformat == FASTQ || CSFASTQ)){
		o.switch2Fasta = true;
		
		if(o.fformat == FASTQ)   o.fformat = FASTA;
		if(o.fformat == CSFASTQ) o.fformat = CSFASTA;
	}
	
	
	getOptionValue(o.readsFile, parser, "reads");
	cout << "Reads file:            " << o.readsFile << endl;
	o.runType = SINGLE;
	
	if(isSet(parser, "reads2")){
		getOptionValue(o.readsFile2, parser, "reads2");
		cout << "Reads file 2:          " << o.readsFile2 << "   (paired run)" << endl;
		o.runType = PAIRED;
	}
	
	
	// barcode and adapter file options
	
	if(isSet(parser, "barcodes")){
		
		getOptionValue(o.barcodeFile, parser, "barcodes");
		cout << "Barcode file:          " << o.barcodeFile << endl;
		
		if(isSet(parser, "barcode-reads")){
			getOptionValue(o.barReadsFile, parser, "barcode-reads");
			cout << "Barcode reads file:    " << o.barReadsFile << endl;
			
			o.barDetect = BARCODE_READ;
		}
		else o.barDetect = WITHIN_READ_REMOVAL;
		
		if(o.runType == SINGLE)      o.runType = SINGLE_BARCODED;
		else if(o.runType == PAIRED) o.runType = PAIRED_BARCODED;
		
		if(o.barDetect == WITHIN_READ_REMOVAL && isSet(parser, "barcode-keep")){
			o.barDetect = WITHIN_READ;
		}
	}
	
	if(isSet(parser, "adapters")){
		getOptionValue(o.adapterFile, parser, "adapters");
		cout << "Adapter file:          " << o.adapterFile << endl;
		o.adapRm = NORMAL;
		o.useAdapterFile = true;
	}
	else if(isSet(parser, "adapter-seq")){
		getOptionValue(o.adapterSeq, parser, "adapter-seq");
		o.adapRm = NORMAL;
	}
	cout << endl;
	
	
	// filtering and trimming options
	
	getOptionValue(o.nThreads, parser, "threads");
	cout << "threads:               " << o.nThreads << endl;
	
	getOptionValue(o.maxUncalled, parser, "max-uncalled");
	cout << "max-uncalled:          " << o.maxUncalled << endl;
	
	if(isSet(parser, "pre-trim-left")){
		getOptionValue(o.cutLen_begin, parser, "pre-trim-left");
		cout << "pre-trim-left:         " << o.cutLen_begin << endl;
	}
	
	if(isSet(parser, "pre-trim-right")){
		getOptionValue(o.cutLen_end, parser, "pre-trim-right");
		cout << "pre-trim-right:        " << o.cutLen_end << endl;
	}
	
	if(isSet(parser, "pre-trim-phred")){
		getOptionValue(o.phred_preQual, parser, "pre-trim-phred");
		cout << "pre-trim-phred:        " << o.phred_preQual << endl;
	}
	
	
	if(isSet(parser, "post-trim-length")){
		getOptionValue(o.cutLen_read, parser, "post-trim-length");
		cout << "post-trim-length:      " << o.cutLen_read << endl;
	}
	
	getOptionValue(o.min_readLen, parser, "min-readlength");
	cout << "min-readlength:        " << o.min_readLen << endl;
	if(o.isColorSpace) o.min_readLen++;
	
	
	// logging and tagging options
	
	if(isSet(parser, "log-level")){
		getOptionValue(o.log_level, parser, "log-level");
		
		     if(o.log_level == "ALL") o.logLevel = ALL;
		else if(o.log_level == "TAB") o.logLevel = TAB;
		else if(o.log_level == "MOD") o.logLevel = MOD;
		else if(o.log_level != "NONE"){
			cerr << "Specified log-level is unknown!" << "\n" << endl;
			exit(1);
		}
	}
	
	if(isSet(parser, "single-reads")) o.writeSingleReads = true;
	if(isSet(parser, "length-dist"))  o.writeLengthDist  = true;
	if(isSet(parser, "number-tags"))  o.useNumberTag     = true;
	if(isSet(parser, "removal-tags")) o.useRemovalTag    = true;
	if(isSet(parser, "random-tags"))  o.randTag          = true;
	
	
	// barcode options
	
	if(o.barDetect != BOFF){
		
		getOptionValue(o.b_trim_end, parser, "barcode-trim-end");
		if     (o.b_trim_end == "LEFT")        o.b_end = LEFT;
		else if(o.b_trim_end == "RIGHT")       o.b_end = RIGHT;
		else if(o.b_trim_end == "ANY")         o.b_end = ANY;
		else if(o.b_trim_end == "LEFT_TAIL")   o.b_end = LEFT_TAIL;
		else if(o.b_trim_end == "RIGHT_TAIL")  o.b_end = RIGHT_TAIL;
		else{
			cerr << "Specified barcode trim-end is unknown!" << "\n" << endl;
			exit(1);
		}
		cout << "barcode-trim-end:      " << o.b_trim_end << endl;
		
		
		if(isSet(parser, "barcode-tail-length")){
			getOptionValue(o.b_tail_len, parser, "barcode-tail-length");
			cout << "barcode-tail-length:   " << o.b_tail_len << endl;
		}
		
		if(isSet(parser, "barcode-min-overlap")){
			getOptionValue(o.b_min_overlap, parser, "barcode-min-overlap");
			cout << "barcode-min-overlap:   " << o.b_min_overlap << endl;
		}
		
		getOptionValue(o.b_threshold, parser, "barcode-threshold");
		cout << "barcode-threshold:     " << o.b_threshold << endl;
		
		getOptionValue(o.b_match, parser, "barcode-match");
		getOptionValue(o.b_mismatch, parser, "barcode-mismatch");
		getOptionValue(o.b_gapCost, parser, "barcode-gap");
		
		cout << "barcode-match:        ";
		if(o.b_match >= 0) cout << " ";
		cout << o.b_match << endl;
		
		cout << "barcode-mismatch:     ";
		if(o.b_mismatch >= 0) cout << " ";
		cout << o.b_mismatch << endl;
		
		cout << "barcode-gap:          ";
		if(o.b_gapCost >= 0) cout << " ";
		cout << o.b_gapCost << "\n" << endl;
	}
	
	
	// adapter options
	
	if(o.adapRm != AOFF){
		
		getOptionValue(o.a_trim_end, parser, "adapter-trim-end");
		if     (o.a_trim_end == "LEFT")        o.end = LEFT;
		else if(o.a_trim_end == "RIGHT")       o.end = RIGHT;
		else if(o.a_trim_end == "ANY")         o.end = ANY;
		else if(o.a_trim_end == "LEFT_TAIL")   o.end = LEFT_TAIL;
		else if(o.a_trim_end == "RIGHT_TAIL")  o.end = RIGHT_TAIL;
		else {
			cerr << "Specified adapter trim-end is unknown!" << "\n" << endl;
			exit(1);
		}
		cout << "adapter-trim-end:      " << o.a_trim_end << endl;
		
		
		if(isSet(parser, "adapter-tail-length")){
			getOptionValue(o.a_tail_len, parser, "adapter-tail-length");
			cout << "adapter-tail-length:   " << o.a_tail_len << endl;
		}
		
		getOptionValue(o.a_min_overlap, parser, "adapter-min-overlap");
		cout << "adapter-min-overlap:   " << o.a_min_overlap << endl;
		
		getOptionValue(o.a_threshold, parser, "adapter-threshold");
		cout << "adapter-threshold:     " << o.a_threshold << endl;
		
		
		getOptionValue(o.match, parser, "adapter-match");
		getOptionValue(o.mismatch, parser, "adapter-mismatch");
		getOptionValue(o.gapCost, parser, "adapter-gap");
		
		cout << "adapter-match:        ";
		if(o.match >= 0) cout << " ";
		cout << o.match << endl;
		
		cout << "adapter-mismatch:     ";
		if(o.mismatch >= 0) cout << " ";
		cout << o.mismatch << endl;
		
		cout << "adapter-gap:          ";
		if(o.gapCost >= 0) cout << " ";
		cout << o.gapCost << "\n" << endl;
	}
	
	
	// option compatibility tests
	
	if(o.cutLen_read != 0 && o.cutLen_read < o.min_readLen){
		o.cutLen_read = 0;
		cerr << "\nOption post-trim-length omitted, as it is shorter than min read length.\n" << endl;
	}
	
}


#endif /* FLEXBAR_OPTIONS_H_ */
