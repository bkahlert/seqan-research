
#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>
 

#include <seqan/seq_io.h>

#include "seqc.h"
#include "read_stats.h"
#include "input.h"
#include "output.h"



using namespace seqan;
using std::cout;
using std::endl;
using std::ostream;


void AppOptions::out(ostream & os){
    os << "active configuration" << endl
       << "========================================"<< endl
       << "VERBOSITY\t" << verbosity << endl
       << "inputFile\t" << inputFile << endl
       << "max-reads\t" << maxReads << endl
       << "outdir\t" << outdir << endl
       << "manual offset\t" << scoreOffset << endl; 
       ;
}


seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("seqc");
    // Set short description, version, and date.
    setShortDescription(parser, "Quality Control app for NGS reads - implemented using seqan lib.");
    setVersion(parser, "0.1");
    setDate(parser, "May 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application description.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "inputFile"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
    
    // custom options
    addOption(parser, seqan::ArgParseOption(
                "r","max-reads",
                "Maximum number of records to analyse. Report on what we found so far, as if the input ended after given number of reads",  
                seqan::ArgParseArgument::INTEGER, "INT")); 

    addOption(parser, seqan::ArgParseOption(
                "o", "outdir",
                "directory path that will be created to contain report data, defaults to <inputFile>_out at inputFile location",
                seqan::ArgParseArgument::STRING, "STRING"
                ));
    addOption(parser, seqan::ArgParseOption(
                "s","score-offset",
                "force the given number as score offset. Don't guess from input data.",  
                seqan::ArgParseArgument::INTEGER, "INT")); 

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


    seqan::getArgumentValue(options.inputFile, parser, 0);
    
    /**
     * set valid outdir
     */
    seqan::getOptionValue(options.outdir, parser, "outdir");
    if(seqan::empty(options.outdir)){
        options.outdir = "out_";
        append(options.outdir, basename(seqan::toCString(options.inputFile)));
        // @TODO strip filename extension
        
    }
    
    seqan::getOptionValue(options.maxReads,parser, "max-reads");
    seqan::getOptionValue(options.scoreOffset,parser, "score-offset");
    
   

    return seqan::ArgumentParser::PARSE_OK;
}

